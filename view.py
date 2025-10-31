"""
DICOMシリーズ（CT/MRI）ビューア

機能:
- ディレクトリ/ファイルからDICOMシリーズ（マルチフレーム含む）読込
- スライダーでスライス番号、Window Center(レベル)/Window Width(幅)の調整
- 断面(Axial/Coronal/Sagittal)の切替

依存: pydicom, PyQt5, numpy
"""

import os
import math
from typing import List, Tuple, Optional

import numpy as np
import pydicom
from pydicom.dataset import FileDataset
from PyQt5 import QtWidgets, QtGui, QtCore


def _safe_get_first(tag_value):
	"""DICOMタグの値が配列/マルチ値でも先頭を返すヘルパー。"""
	if tag_value is None:
		return None
	if isinstance(tag_value, (list, tuple)):
		return tag_value[0] if len(tag_value) > 0 else None
	return tag_value


def windowing(image: np.ndarray, center: float, width: float, invert: bool = False) -> np.ndarray:
	"""Window Level/Widthを適用して8bitグレースケール画像に変換。

	Args:
		image: 2D numpy array (float or int)
		center: window center (level)
		width: window width
		invert: MONOCHROME1の場合にTrue
	Returns:
		uint8 2D array in [0,255]
	"""
	if width is None or width <= 1e-6:
		# 幅がゼロに近い場合は自動設定
		vmin, vmax = np.nanmin(image), np.nanmax(image)
	else:
		vmin = center - width / 2.0
		vmax = center + width / 2.0
	img = np.clip(image.astype(np.float32), vmin, vmax)
	if vmax - vmin > 0:
		img = (img - vmin) / (vmax - vmin)
	else:
		img = np.zeros_like(img, dtype=np.float32)
	img8 = (img * 255.0 + 0.5).astype(np.uint8)
	if invert:
		img8 = 255 - img8
	return img8


def to_qimage_gray(img8: np.ndarray) -> QtGui.QImage:
	"""uint8 2D配列をQImage(Grayscale8)へ。"""
	h, w = img8.shape
	# bytesPerLine は幅に等しい（1byte/pixel）
	qimg = QtGui.QImage(img8.data, w, h, w, QtGui.QImage.Format_Grayscale8)
	# バッファがガベージコレクションで解放されないように参照を保持
	qimg.ndarray = img8  # type: ignore
	return qimg


class DicomVolume:
	"""DICOMシリーズを3Dボリュームにまとめたデータ構造。

	Attributes:
		volume: 3D numpy array (Z, Y, X)
		spacing: (z, y, x) mm
		photometric: 'MONOCHROME1' or 'MONOCHROME2'
		default_window: (center, width) or None
	"""

	def __init__(self, volume: np.ndarray, spacing: Tuple[float, float, float],
				 photometric: str = "MONOCHROME2",
				 default_window: Optional[Tuple[float, float]] = None) -> None:
		self.volume = volume.astype(np.float32)
		self.spacing = spacing
		self.photometric = photometric
		self.default_window = default_window


def _apply_rescale(arr: np.ndarray, ds: FileDataset) -> np.ndarray:
	"""Rescale Slope/Intercept を適用。"""
	slope = float(getattr(ds, 'RescaleSlope', 1) or 1)
	intercept = float(getattr(ds, 'RescaleIntercept', 0) or 0)
	return arr.astype(np.float32) * slope + intercept


def _sort_key_for_slice(ds: FileDataset) -> float:
	# ImagePositionPatientがあればzでソート、なければInstanceNumber、それもなければ0
	ipp = getattr(ds, 'ImagePositionPatient', None)
	if ipp is not None and len(ipp) >= 3:
		try:
			return float(ipp[2])
		except Exception:
			pass
	inst = getattr(ds, 'InstanceNumber', None)
	if inst is not None:
		try:
			return float(inst)
		except Exception:
			pass
	# 最後の手段
	return 0.0


def load_dicom_series(path: str) -> DicomVolume:
	"""パス（ファイル or ディレクトリ）からDICOMシリーズを読み込み、3Dボリュームにする。

	- ディレクトリ: 同一SeriesInstanceUIDごとにグループ化し、枚数最大のシリーズを採用
	- ファイル: マルチフレームは(Frames, Rows, Cols)として扱う。単一フレームは1枚ボリューム
	"""
	if os.path.isdir(path):
		dcms: List[FileDataset] = []
		for root, _, files in os.walk(path):
			for fn in files:
				if fn.lower().endswith(('.dcm', '.dicom')) or '.' not in fn:
					fp = os.path.join(root, fn)
					try:
						ds = pydicom.dcmread(fp, force=True, stop_before_pixels=False)
						# ピクセルデータを持つもののみ
						if 'PixelData' in ds:
							dcms.append(ds)
					except Exception:
						pass
		if not dcms:
			raise RuntimeError('DICOMファイルが見つかりませんでした。')

		# SeriesInstanceUIDでグループ化
		series_map = {}
		for ds in dcms:
			sid = getattr(ds, 'SeriesInstanceUID', 'unknown-series')
			series_map.setdefault(sid, []).append(ds)

		# 最大枚数のシリーズを選択
		sid, series = max(series_map.items(), key=lambda kv: len(kv[1]))

		# スライス順に並べ替え
		series.sort(key=_sort_key_for_slice)

		# spacing
		# In-plane: PixelSpacing (row, col) = (y, x)
		ps = getattr(series[0], 'PixelSpacing', [1.0, 1.0])
		try:
			yx = (float(ps[0]), float(ps[1]))
		except Exception:
			yx = (1.0, 1.0)

		# z方向の間隔: SpacingBetweenSlices or SliceThickness or 差分
		dz = None
		try:
			dz = float(getattr(series[0], 'SpacingBetweenSlices'))
		except Exception:
			try:
				dz = float(getattr(series[0], 'SliceThickness'))
			except Exception:
				pass
		if dz is None and len(series) >= 2:
			z_positions = []
			for ds in series:
				ipp = getattr(ds, 'ImagePositionPatient', None)
				if ipp is not None and len(ipp) >= 3:
					try:
						z_positions.append(float(ipp[2]))
					except Exception:
						pass
			if len(z_positions) >= 2:
				diffs = np.diff(sorted(z_positions))
				if diffs.size > 0:
					dz = float(np.median(np.abs(diffs)))
		if dz is None:
			dz = 1.0

		photometric = getattr(series[0], 'PhotometricInterpretation', 'MONOCHROME2')

		# default WL/WW
		wc = _safe_get_first(getattr(series[0], 'WindowCenter', None))
		ww = _safe_get_first(getattr(series[0], 'WindowWidth', None))
		default_window = None
		try:
			if wc is not None and ww is not None:
				default_window = (float(wc), float(ww))
		except Exception:
			default_window = None

		# 3D volume生成
		slices = []
		for ds in series:
			arr = ds.pixel_array
			# 2D想定
			if arr.ndim == 3 and arr.shape[0] > 1:
				# まれにマルチフレームが混ざるのを防ぐため先頭のみ
				arr = arr[0]
			arr = _apply_rescale(arr, ds)
			slices.append(arr)
		vol = np.stack(slices, axis=0)  # (Z, Y, X)

		return DicomVolume(vol, spacing=(float(dz), float(yx[0]), float(yx[1])),
						   photometric=photometric, default_window=default_window)

	else:
		# 単一ファイル
		ds = pydicom.dcmread(path, force=True, stop_before_pixels=False)
		if 'PixelData' not in ds:
			raise RuntimeError('選択したファイルにPixelDataが含まれていません。')
		arr = ds.pixel_array
		photometric = getattr(ds, 'PhotometricInterpretation', 'MONOCHROME2')

		wc = _safe_get_first(getattr(ds, 'WindowCenter', None))
		ww = _safe_get_first(getattr(ds, 'WindowWidth', None))
		default_window = None
		try:
			if wc is not None and ww is not None:
				default_window = (float(wc), float(ww))
		except Exception:
			default_window = None

		# spacing
		ps = getattr(ds, 'PixelSpacing', [1.0, 1.0])
		try:
			yx = (float(ps[0]), float(ps[1]))
		except Exception:
			yx = (1.0, 1.0)
		dz = float(getattr(ds, 'SpacingBetweenSlices', getattr(ds, 'SliceThickness', 1.0)) or 1.0)

		if arr.ndim == 3:  # Multi-frame (Frames, Rows, Cols)
			arr = _apply_rescale(arr, ds)
			vol = arr.astype(np.float32)
		elif arr.ndim == 2:
			vol = _apply_rescale(arr, ds)[None, ...]  # (1, Y, X)
		else:
			raise RuntimeError(f'未対応のpixel_array形状: {arr.shape}')

		return DicomVolume(vol, spacing=(float(dz), float(yx[0]), float(yx[1])),
						   photometric=photometric, default_window=default_window)


class ImageView(QtWidgets.QLabel):
	"""画像表示用ラベル。リサイズ時に自動でフィット表示。"""

	def __init__(self):
		super().__init__()
		self.setAlignment(QtCore.Qt.AlignCenter)
		self._pixmap: Optional[QtGui.QPixmap] = None
		# 右下オーバーレイ用（プレビュー等）
		self._overlay_widget: Optional[QtWidgets.QWidget] = None
		self._overlay_margin = 8

	def set_image(self, qimg: QtGui.QImage):
		self._pixmap = QtGui.QPixmap.fromImage(qimg)
		self._update_scaled_pixmap()

	def clear_image(self):
		self._pixmap = None
		self.clear()

	def set_overlay_widget(self, w: Optional[QtWidgets.QWidget]):
		"""このビューの右下に重ねる子ウィジェットを設定。"""
		if self._overlay_widget is not None:
			self._overlay_widget.setParent(None)
		self._overlay_widget = w
		if w is not None:
			w.setParent(self)
			w.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents, True)
			w.show()
			self._position_overlay()

	def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
		super().resizeEvent(event)
		self._update_scaled_pixmap()
		self._position_overlay()

	def _update_scaled_pixmap(self):
		if self._pixmap is None:
			return
		scaled = self._pixmap.scaled(self.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
		self.setPixmap(scaled)

	def _position_overlay(self):
		if self._overlay_widget is None:
			return
		m = self._overlay_margin
		over = self._overlay_widget
		# 右下に配置
		sz = over.size()
		x = max(0, self.width() - sz.width() - m)
		y = max(0, self.height() - sz.height() - m)
		over.move(x, y)


class MainWindow(QtWidgets.QMainWindow):
	def __init__(self):
		super().__init__()
		self.setWindowTitle('DICOM Viewer (CT/MRI)')
		self.resize(1100, 800)

		# 状態
		self.volume: Optional[DicomVolume] = None
		self.orientation = 'Axial'  # 'Axial'|'Coronal'|'Sagittal'
		self.slice_index = 0
		self.window_center: Optional[float] = None
		self.window_width: Optional[float] = None

		# UI
		central = QtWidgets.QWidget()
		self.setCentralWidget(central)
		hbox = QtWidgets.QHBoxLayout(central)

		self.image_view = ImageView()
		hbox.addWidget(self.image_view, stretch=1)

		# 右下Axialプレビュー（ウィンドウに重ねる）
		self.axial_preview = QtWidgets.QLabel(self)
		self.axial_preview.setFixedSize(220, 220)
		self.axial_preview.setFrameShape(QtWidgets.QFrame.Box)
		self.axial_preview.setLineWidth(1)
		self.axial_preview.setStyleSheet('background-color: rgba(0,0,0,0);')
		self.axial_preview.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents, True)
		self.axial_preview.setAlignment(QtCore.Qt.AlignCenter)
		self.axial_preview.setToolTip('Axial preview (mid-slice)')
		self.axial_preview.raise_()
		self._position_axial_preview()

		# コントロールパネル
		panel = QtWidgets.QFrame()
		panel.setFrameShape(QtWidgets.QFrame.StyledPanel)
		panel.setMaximumWidth(360)
		v = QtWidgets.QVBoxLayout(panel)

		# 開くボタン
		open_layout = QtWidgets.QHBoxLayout()
		btn_open_dir = QtWidgets.QPushButton('フォルダを開く')
		btn_open_file = QtWidgets.QPushButton('ファイルを開く')
		open_layout.addWidget(btn_open_dir)
		open_layout.addWidget(btn_open_file)
		v.addLayout(open_layout)

		# 断面切替
		form = QtWidgets.QFormLayout()
		self.combo_orientation = QtWidgets.QComboBox()
		self.combo_orientation.addItems(['Axial', 'Coronal', 'Sagittal'])
		form.addRow('断面', self.combo_orientation)

		# スライススライダー
		self.slider_slice = QtWidgets.QSlider(QtCore.Qt.Horizontal)
		self.slider_slice.setMinimum(0)
		self.slider_slice.setMaximum(0)
		self.slider_slice.setValue(0)
		self.slider_slice.setEnabled(False)
		self.lbl_slice = QtWidgets.QLabel('0 / 0')
		form.addRow('スライス', self.slider_slice)
		form.addRow('位置', self.lbl_slice)

		# Window level/width
		self.slider_wc = QtWidgets.QSlider(QtCore.Qt.Horizontal)
		self.slider_ww = QtWidgets.QSlider(QtCore.Qt.Horizontal)
		for s in (self.slider_wc, self.slider_ww):
			s.setMinimum(-4096)
			s.setMaximum(4096)
			s.setEnabled(False)
		self.lbl_wc = QtWidgets.QLabel('WC: -')
		self.lbl_ww = QtWidgets.QLabel('WW: -')
		form.addRow('Window Center', self.slider_wc)
		form.addRow(self.lbl_wc)
		form.addRow('Window Width', self.slider_ww)
		form.addRow(self.lbl_ww)

		v.addLayout(form)
		v.addStretch(1)

		hbox.addWidget(panel, stretch=0)

		# ステータスバー
		self.status = QtWidgets.QStatusBar()
		self.setStatusBar(self.status)
		self._update_status()

		# シグナル
		btn_open_dir.clicked.connect(self._on_open_dir)
		btn_open_file.clicked.connect(self._on_open_file)
		self.combo_orientation.currentTextChanged.connect(self._on_orientation_changed)
		self.slider_slice.valueChanged.connect(self._on_slice_changed)
		self.slider_wc.valueChanged.connect(self._on_wc_changed)
		self.slider_ww.valueChanged.connect(self._on_ww_changed)

	def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
		"""ウィンドウのリサイズ時にプレビューを右下へ再配置。"""
		super().resizeEvent(event)
		self._position_axial_preview()

	# --- UI handlers ---
	def _on_open_dir(self):
		path = QtWidgets.QFileDialog.getExistingDirectory(self, 'DICOMフォルダを選択')
		if not path:
			return
		self._load_path(path)

	def _on_open_file(self):
		path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'DICOMファイルを選択', filter='DICOM (*.dcm *.dicom);;All (*)')
		if not path:
			return
		self._load_path(path)

	def _on_orientation_changed(self, text: str):
		self.orientation = text
		self._reset_slice_slider()
		self._update_image()

	def _on_slice_changed(self, value: int):
		self.slice_index = value
		self._update_image()

	def _on_wc_changed(self, value: int):
		self.window_center = float(value)
		self.lbl_wc.setText(f'WC: {self.window_center:.1f}')
		self._update_image()

	def _on_ww_changed(self, value: int):
		# WWは正にする（ゼロや負だと見えにくい）
		self.window_width = float(max(1, abs(value)))
		self.lbl_ww.setText(f'WW: {self.window_width:.1f}')
		self._update_image()

	# --- helpers ---
	def _load_path(self, path: str):
		try:
			vol = load_dicom_series(path)
		except Exception as e:
			QtWidgets.QMessageBox.critical(self, '読み込みエラー', str(e))
			return

		self.volume = vol
		# デフォルトのウィンドウ
		if vol.default_window is not None:
			self.window_center, self.window_width = vol.default_window
		else:
			vmin, vmax = float(np.nanmin(vol.volume)), float(np.nanmax(vol.volume))
			self.window_center = (vmin + vmax) / 2.0
			self.window_width = max(1.0, vmax - vmin)

		# スライダー有効化
		self.slider_wc.setEnabled(True)
		self.slider_ww.setEnabled(True)
		self.slider_wc.setValue(int(self.window_center))
		self.slider_ww.setValue(int(self.window_width))

		self._reset_slice_slider()
		self._update_image()
		self._update_axial_preview()
		self._update_status()

	def _reset_slice_slider(self):
		if self.volume is None:
			self.slider_slice.setEnabled(False)
			self.slider_slice.setMinimum(0)
			self.slider_slice.setMaximum(0)
			self.slider_slice.setValue(0)
			self.lbl_slice.setText('0 / 0')
			return
		z, y, x = self.volume.volume.shape
		if self.orientation == 'Axial':
			max_idx = z - 1
		elif self.orientation == 'Coronal':
			max_idx = y - 1
		else:
			max_idx = x - 1
		max_idx = max(0, int(max_idx))
		self.slider_slice.setEnabled(True)
		self.slider_slice.setMinimum(0)
		self.slider_slice.setMaximum(max_idx)
		# 現在位置を近い値に
		self.slice_index = min(self.slice_index, max_idx)
		self.slider_slice.setValue(self.slice_index)
		self.lbl_slice.setText(f'{self.slice_index + 1} / {max_idx + 1}')

	def _get_current_slice2d(self) -> Optional[np.ndarray]:
		if self.volume is None:
			return None
		vol = self.volume.volume
		idx = int(self.slice_index)
		# 断面に応じて2Dスライスを抽出
		if self.orientation == 'Axial':
			idx = np.clip(idx, 0, vol.shape[0] - 1)
			sl = vol[idx, :, :]
		elif self.orientation == 'Coronal':
			idx = np.clip(idx, 0, vol.shape[1] - 1)
			sl = vol[:, idx, :]
		else:  # Sagittal
			idx = np.clip(idx, 0, vol.shape[2] - 1)
			sl = vol[:, :, idx]
		return sl

	def _update_image(self):
		sl = self._get_current_slice2d()
		if sl is None:
			self.image_view.clear_image()
			# プレビューもクリア
			self.axial_preview.clear()
			return
		wc = self.window_center if self.window_center is not None else float(np.mean(sl))
		ww = self.window_width if self.window_width is not None else float(np.ptp(sl))
		invert = (str(self.volume.photometric).upper() == 'MONOCHROME1') if self.volume else False
		img8 = windowing(sl, wc, ww, invert=invert)
		qimg = to_qimage_gray(img8)
		# 物理縦横比を維持するため、行方向/列方向の間隔(dy, dx)を用いて非等方スケーリング
		if self.volume is not None:
			rows, cols = img8.shape
			# 断面に応じた平面内の物理ピクセル間隔を取得
			dy, dx = self._get_plane_spacing()
			# 不正な値を避ける
			if dx is None or dx <= 0:
				dx = dy if (dy is not None and dy > 0) else 1.0
			if dy is None or dy <= 0:
				dy = 1.0
			# 各ピクセルの縦横比をdy/dxに合わせるため、画像の高さを倍率(dy/dx)でスケーリング
			target_w = int(cols)
			target_h = max(1, int(round(rows * (dy / dx))))
			# メモリ過負荷を避ける簡易キャップ
			cap = 4096
			m = max(target_w, target_h)
			if m > cap:
				scale = cap / float(m)
				target_w = max(1, int(round(target_w * scale)))
				target_h = max(1, int(round(target_h * scale)))
			qimg = qimg.scaled(target_w, target_h, QtCore.Qt.IgnoreAspectRatio, QtCore.Qt.SmoothTransformation)
		self.image_view.set_image(qimg)
		# ラベル更新
		if self.slider_slice.isEnabled():
			max_idx = self.slider_slice.maximum()
			self.lbl_slice.setText(f'{self.slice_index + 1} / {max_idx + 1}')
		self._update_status()
		# プレビュー更新
		self._update_axial_preview()

	def _update_axial_preview(self):
		"""右下のAxialプレビューを更新。Sagittal/Coronal時は位置ラインを重畳。"""
		if self.volume is None or self.axial_preview is None:
			return
		vol = self.volume.volume
		z = vol.shape[0]
		if z <= 0:
			self.axial_preview.clear()
			return
		mid_idx = max(0, min(z - 1, z // 2))
		sl = vol[mid_idx, :, :]
		# WC/WW（未設定時はスライスから推定）
		wc = self.window_center if self.window_center is not None else float(np.mean(sl))
		ww = self.window_width if self.window_width is not None else float(np.ptp(sl))
		invert = (str(self.volume.photometric).upper() == 'MONOCHROME1')
		img8 = windowing(sl, wc, ww, invert=invert)
		qimg = to_qimage_gray(img8)

		# Axial平面の物理比で縦横比補正
		rows, cols = img8.shape
		_, sy, sx = self.volume.spacing
		target_w = int(cols)
		target_h = max(1, int(round(rows * (sy / sx if sx > 0 else 1.0))))
		qimg_phys = qimg.scaled(target_w, target_h, QtCore.Qt.IgnoreAspectRatio, QtCore.Qt.SmoothTransformation)

		# プレビューサイズにレターボックスでフィット
		pw = self.axial_preview.width()
		ph = self.axial_preview.height()
		canvas = QtGui.QImage(pw, ph, QtGui.QImage.Format_ARGB32)
		canvas.fill(QtGui.QColor(0, 0, 0, 0))
		painter = QtGui.QPainter(canvas)
		try:
			# 背面を半透明黒に
			painter.fillRect(0, 0, pw, ph, QtGui.QColor(0, 0, 0, 120))
			# スケーリング計算
			src_w = qimg_phys.width()
			src_h = qimg_phys.height()
			scale = min(pw / src_w, ph / src_h) if src_w > 0 and src_h > 0 else 1.0
			dw = int(round(src_w * scale))
			dh = int(round(src_h * scale))
			offx = (pw - dw) // 2
			offy = (ph - dh) // 2
			# 画像描画
			target_rect = QtCore.QRect(offx, offy, dw, dh)
			painter.drawImage(target_rect, qimg_phys)

			# 現在の断面に応じてラインを描画
			pen_sag = QtGui.QPen(QtGui.QColor(255, 80, 80), 2)
			pen_cor = QtGui.QPen(QtGui.QColor(80, 255, 80), 2)
			if self.orientation == 'Sagittal':
				# x方向位置
				x_dim = vol.shape[2]
				fx = 0.0 if x_dim <= 1 else float(self.slice_index) / float(x_dim - 1)
				px = offx + int(round(fx * dw))
				painter.setPen(pen_sag)
				painter.drawLine(px, offy, px, offy + dh)
			elif self.orientation == 'Coronal':
				# y方向位置
				y_dim = vol.shape[1]
				fy = 0.0 if y_dim <= 1 else float(self.slice_index) / float(y_dim - 1)
				py = offy + int(round(fy * dh))
				painter.setPen(pen_cor)
				painter.drawLine(offx, py, offx + dw, py)
		finally:
			painter.end()

		self.axial_preview.setPixmap(QtGui.QPixmap.fromImage(canvas))

	def _get_plane_spacing(self) -> Tuple[float, float]:
		"""現在の断面に対応する2D平面の(行方向dy, 列方向dx)の物理間隔(mm)を返す。

		Axial:   rows=Y(cols= X) -> (sy, sx)
		Coronal: rows=Z(cols= X) -> (sz, sx)
		Sagittal:rows=Z(cols= Y) -> (sz, sy)
		"""
		if self.volume is None:
			return (1.0, 1.0)
		sz, sy, sx = self.volume.spacing
		if self.orientation == 'Axial':
			return (sy, sx)
		elif self.orientation == 'Coronal':
			return (sz, sx)
		else:  # Sagittal
			return (sz, sy)

	def _update_status(self):
		if self.volume is None:
			self.status.showMessage('準備完了 - ファイル/フォルダからDICOMを開いてください')
			return
		z, y, x = self.volume.volume.shape
		sz, sy, sx = self.volume.spacing
		msg = f'Volume: {z}x{y}x{x} | Spacing(mm): {sz:.3f}/{sy:.3f}/{sx:.3f} | Ori: {self.orientation} | Slice: {self.slice_index + 1} '
		if self.window_center is not None and self.window_width is not None:
			msg += f'| WC/WW: {self.window_center:.1f}/{self.window_width:.1f}'
		self.status.showMessage(msg)

	def _position_axial_preview(self):
		"""QMainWindow の内容領域右下へプレビューを配置。"""
		if getattr(self, 'axial_preview', None) is None:
			return
		m = 8
		rect = self.contentsRect()
		top_left = rect.topLeft()
		w = self.axial_preview.width()
		h = self.axial_preview.height()
		x = top_left.x() + rect.width() - w - m
		y = top_left.y() + rect.height() - h - m
		self.axial_preview.move(max(0, x), max(0, y))


def main():
	import sys
	app = QtWidgets.QApplication(sys.argv)
	w = MainWindow()
	w.show()
	sys.exit(app.exec_())


if __name__ == '__main__':
	main()

