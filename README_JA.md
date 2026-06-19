# ORCA Hessian → VEDA .fmu 変換ツール

現在のリリース: **v1.1.0**

## 作者

- 山本 薫（岡山理科大学）

このリポジトリは、**ORCA** の振動解析結果を **VEDA** で扱うための、外部依存の少ない小さなツール集です。

- `orca_hess_to_veda_fmu.py` / `.pyw`  
  ORCA の `*.out` と対応する `*.hess` から、VEDA が読める `*.fmu`
  （Gaussian の fchk 互換“最小”テキスト）を生成します。既定では補助的な `*.mpo` 結合情報ファイルも生成します。

- `adjust_mpo.py` / `.pyw`  
  `.mpo` の **生成／再構成（修復）** に特化したツールです。すでに `.fmu`/`.xyz` がある場合や、
  金属錯体などで結合定義の手動調整が必要な場合に使用します。

English documentation: see [`README.md`](README.md).

## 動作環境

- Python 3.8 以上（推奨）
- メイン変換GUIは標準ライブラリのみで動作します（GUI は `tkinter` を使用）
- 任意: `tkinterdnd2`（ウィンドウ内ドラッグ&ドロップ用）
- 任意: `numpy`（Hessian/周波数の健全性チェック用。変換自体は `numpy` なしでも可能です）

## 使い方（GUI）

1. Windows では `orca_hess_to_veda_fmu.pyw` をダブルクリック、または次を実行します。
   ```bash
   python orca_hess_to_veda_fmu.py
   ```
2. ORCA の `molecule.out` を選択します。
3. 対応する ORCA の `molecule.hess` を選択します。
4. 出力先 `molecule.fmu` を確認します。
5. **Convert** をクリックします。
6. 出力先と同じフォルダに以下が出力されます。
   - `molecule.fmu`
   - 既定では `molecule.mpo`、既存の `.mpo` を保持する場合は `molecule_gen.mpo`

## 変換アルゴリズムの概要

1. ORCA `.out` の最後の **VIBRATIONAL FREQUENCIES** セクションを検出します。
2. その周波数セクションに最も近いCartesian座標ブロックを選択します。
3. 選択された座標ブロックの近傍から電荷と多重度を読み取ります。
4. `.hess` の `$atoms` ブロックが利用できる場合は、その座標を優先します。これはHessianに対応する直接的な座標情報であるためです。
5. `.hess` の `$hessian` ブロックを読み取り、完全なCartesian Hessian行列を復元します。
6. 必要に応じて、質量重み付きHessianの解除やHessian単位のHartree/Bohr²への変換を行います。
7. 原子番号、座標、Cartesian force constantsを含む最小fchk互換 `.fmu` を出力します。
8. 必要に応じて、保守的な距離判定に基づく補助 `.mpo` を生成します。

## `.mpo`（結合情報）の扱い

VEDAは結合を自動推定できますが、金属錯体では過剰結合が起こり、環構造が過剰に生成されることがあります。そのため本ツールが生成する `.mpo` は保守的なヒューリスティックを使います。

- 非金属間結合は共有結合半径と距離マージンから推定します。
- 金属–配位子結合は典型的なドナー原子（`N`, `O`, `S`, `P`, ハロゲン）に制限します。
- 金属の配位は最大配位数と距離カットオフで制限します。
- 既存の `.mpo` は確認なしに上書きしません。

難しい系では、`adjust_mpo.py` / `.pyw` で `.mpo` を再生成または修復してください。

## VEDA 側

VEDA は通常 `.fmu` / `.fmt` を要求します。

- 生成された `molecule.fmu` を指定してください。
- ワークフローによっては `.fmt`（ログ相当）が必要な場合があります。その場合は、ORCA の `molecule.out` をコピーして `molecule.fmt` にリネームし、同じフォルダへ置くことで対応できる場合があります。
- VEDAが結合情報を要求する場合は、生成された `.mpo` を指定し、必要に応じて内容を確認・編集してください。

## ファイル構成

- `orca_hess_to_veda_fmu.py` — コンソール対応GUIエントリポイント
- `orca_hess_to_veda_fmu.pyw` — Windows向け no-console GUIエントリポイント
- `adjust_mpo.py` — コンソール対応 `.mpo` 生成／修復ツール
- `adjust_mpo.pyw` — Windows向け no-console GUIエントリポイント
- `docs/MANUAL_EN.md` — 詳細マニュアル（英語）
- `docs/MANUAL_JA.md` — 詳細マニュアル（日本語）
- `CHANGELOG.md` — 変更履歴
- `CITATION.cff` — 引用メタデータ
- `LICENSE` — MIT License

## 注意点

- 本ツールが作る `.fmu` は Gaussian の `.log` ではありません。VEDAが必要とする **原子情報・座標・Cartesian force constants（Hessian）** を含む最小のfchk互換テキストです。
- `.mpo` は距離ベースの簡易推定です。金属錯体など難しい系では、生成された `.mpo` を確認し、必要に応じて `adjust_mpo.py` を使用してください。
- 周波数チェックは診断用であり、ORCA本体の振動解析を置き換えるものではありません。

## 引用

このプログラムを論文中の結果生成に使用した場合は、計算に使用した **正確なリリース版** を引用してください。

推奨手順:

1. `CITATION.cff` にリリース版、公開日、リポジトリURL、利用可能であればDOIを記載します。
2. GitHubで `v1.1.0` のようなリリースタグを作成します。
3. ZenodoなどDOIを発行するアーカイブにそのリリースを保存します。
4. 論文ではアーカイブ済みリリースのDOIを引用し、計算詳細にはバージョン番号も明記します。

論文中の記載例:

> ORCA Hessian output was converted to VEDA-readable `.fmu` files using ORCA Hessian → VEDA .fmu Converter v1.1.0 (Yamamoto, 2026).

## サポート

- 不具合報告・要望は **GitHub Issues** にお願いします。公開で共有されることで、他のユーザーの再現・解決にも役立ちます。
- 非公開での連絡が必要な場合は、以下へお願いします：  
  `k-yamamoto@ous.ac.jp`

## ライセンス

MIT License（`LICENSE` を参照）。
