# DICOM Viewer

PyQt5 + pydicom によるDICOMシリーズビューア

![使用した際の画像](./画像2.png)

## 機能
- DICOMシリーズ（ディレクトリ）または単一ファイルを読み込み
- スライダーでスライス位置を変更
- Window Center/Width をスライダーで調整
- 断面（Axial / Coronal / Sagittal）切り替え

## 使用方法
### Windowsの場合
- `viewer.exe`を起動

### Linux等の場合(ソースコードから実行する場合)

⚠️注意: ディレクトリ名に日本語が含まれている場合、実行できない可能性があります。

#### 依存パッケージのインストール
```
pip install -r requirements.txt
```

#### 起動方法
```
python view.py
```

起動後、右上の「フォルダを開く」または「ファイルを開く」からDICOMデータを選択してください。

## 使い方
- 断面: Axial/Coronal/Sagittal を切替
- スライス: スライダーでスライス番号を変更
- Window Center/Width: スライダーで明るさ・コントラストを調整

## 注意
- 異なるシリーズが混在するフォルダを選ぶと、枚数が最大のシリーズを自動選択します。
- PhotometricInterpretation が MONOCHROME1 の場合は白黒を反転させます。
- 一部のDICOMではタグ不足により間隔やWL/WWが推定となることがあります。
