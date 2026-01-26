# ORCA Hessian → VEDA FMU 変換ツール

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18377035.svg)](https://doi.org/10.5281/zenodo.18377035)

ORCA の Hessian ファイル（`.hess`）から、VEDA が読み込める `.fmu` ファイル（Gaussian の formatted checkpoint 互換の最小テキスト形式）を生成します。  
必要に応じて、VEDA 用の結合補助ファイル `.mpo` も生成できます。

---

## 主な機能

- ORCA `.hess` → VEDA `.fmu`（最小 fchk互換テキスト）へ変換
- 必要に応じて `.mpo`（結合定義補助）を生成
- Windows向けGUIランチャ（`.pyw`）を同梱
- 英語／日本語マニュアル同梱
- サンプル入出力（`examples/`）同梱

---

## 動作環境

- Python 3.x
- 標準ライブラリのみ（多くの環境で tkinter も標準搭載）

> Linux環境によっては `tkinter` の別途インストールが必要な場合があります。

---

## インストール

### 方法A：GitHubからZIPで取得
1. リポジトリをZIPでダウンロード
2. 展開する
3. 下記の方法で実行

### 方法B：git clone
```bash
git clone https://github.com/fatalfailure/orca_hess_to_veda_fmu.git
cd orca_hess_to_veda_fmu
```

---

## 使い方（GUI）

### 変換ツール（メイン）
- Windows：  
  `orca_hess_to_veda_fmu.pyw` をダブルクリック
- またはPythonから起動：
```bash
python orca_hess_to_veda_fmu.py
```

手順：
1. ORCA の `.hess` を選択
2. （任意）`.out` を選択（charge/multiplicity を取得）
3. （任意）`.mpo` 生成をON
4. **Convert** を押す
5. 同じフォルダに以下が生成されます：
   - `*.fmu`（必ず生成）
   - `*.mpo`（任意）

---

## 使い方（.mpo のみ作る場合）

`.mpo` を作り直したい場合は：

- Windows：`build_mpo.pyw`
- Python実行：
```bash
python build_mpo.py
```

---

## 出力について

### `.fmu`
VEDA 用の「Gaussian formatted checkpoint（fchk）互換の最小テキスト」です。  
以下の情報を含みます：

- 原子数
- 電荷 / 多重度（取得できる場合）
- 原子番号
- Cartesian 座標
- Cartesian 力定数（ヘッセ行列、下三角パック形式）

### `.mpo`（任意）
VEDA内部処理（例：DD2構築）を補助する結合定義ファイルです。  
距離ベースの簡易判定により生成されます。

---

## サンプル

`examples/` フォルダにサンプル入出力と簡単な使用例があります。

---

## 引用（Cite）

本ソフトウェアを学術研究で使用した場合は、以下を引用してください：

山本 薫, *ORCA Hessian to VEDA FMU Converter* (v1.0.1), Zenodo, https://doi.org/10.5281/zenodo.18377035

---

## 作者

- 山本 薫（岡山理科大学）

---

## サポート

- 不具合報告・要望は **GitHub Issues** にお願いします（公開で共有され、他のユーザーの助けにもなります）。
- 非公開での連絡が必要な場合：  
  `k-yamamoto@ous.ac.jp`

---

## ライセンス

MIT License で公開しています。  
詳しくは `LICENSE` を参照してください。
