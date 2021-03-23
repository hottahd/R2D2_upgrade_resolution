# R2D2_upgrade_resolution

R2D2のデータをアップグレードするためのルーチン
- 解像度変更
- 領域変更

をおこなう.

### パラメタ設定

`namelist.in`で以下のパラメタを設定できる.
-　`caseid` : 元の`caseid`
- `caseid_out` : アップグレードしたものを保存する`caseid`
- `c_nd_read` : 読み込む時間ステップ. 文字列. 最後のステップを読み込みたい時は`end_step`