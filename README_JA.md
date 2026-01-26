# ORCA Hessian → VEDA .fmu 変換ツール

このリポジトリは、**ORCA** の振動解析結果を **VEDA** で扱うための、外部依存のない小さなツール集です。

- `orca_hess_to_veda_fmu.py` / `.pyw`  
  ORCA の `*.hess`（必要に応じて `*.out`）から、VEDA が読める `*.fmu`
  （Gaussian の fchk 互換“最小”テキスト）を生成します。オプションで `*.mpo` も生成できます。

- `build_mpo.py` / `.pyw`  
  `.mpo` の **生成／再構成（修復）** に特化したツールです。すでに `.fmu`/`.xyz` がある場合や、
  金属錯体などで結合定義の調整が必要な場合に便利です。

English documentation: see [`README.md`](README.md).

## 動作環境
- Python 3.8 以上（推奨）
- 標準ライブラリのみ（GUI は `tkinter` を使用）

## 使い方

注: Hessianに関する上級者向け設定は **Advanced settings（高度な設定）** に折り畳んであります。通常はデフォルトのままで問題ありません。（GUI）
1. Windows なら `orca_hess_to_veda_fmu.pyw` をダブルクリック  
   あるいは次のように実行：
   ```bash
   python orca_hess_to_veda_fmu.py
   ```
2. `molecule.hess` を選択（必須）
3. （任意）`molecule.out` を選択
4. **Convert** をクリック
5. 入力ファイルと同じフォルダに出力：
   - `molecule.fmu`
   - `molecule.mpo`（任意）

## VEDA 側
VEDA は通常 `.fmu` / `.fmt` を要求します。

- 生成された `molecule.fmu` を指定してください。
- ワークフローによっては `.fmt`（ログ相当）が必要な場合があります。その場合は、
  `molecule.out` をコピーして `molecule.fmt` にリネームして同じフォルダへ置く、などで対応できます。

## ドキュメント
- `docs/MANUAL_EN.md` — 詳細マニュアル（英語）
- `docs/MANUAL_JA.md` — 詳細マニュアル（日本語）

## 注意点
- 本ツールが作る `.fmu` は Gaussian の `.log` ではありません。VEDAが必要とする
  **原子情報・座標・ヘッセ行列**を含む最小の fchk 互換テキストです。
- `.mpo` は距離ベースの簡易推定です。金属錯体など難しい系では `build_mpo.py` による
  生成／修復を推奨します。

## ライセンス
MIT License（`LICENSE` を参照）。
