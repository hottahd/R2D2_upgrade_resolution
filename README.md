# R2D2_upgrade_resolution

R2D2のデータをアップグレードするためのルーチン
- 解像度変更
- 領域変更

をおこなう.

### パラメタ設定

`namelist.in`で以下のパラメタを設定できる.

- `caseid` : 元の`caseid`
- `caseid_out` : アップグレードしたものを保存する`caseid`
- `c_nd_read` : 読み込む時間ステップ. 文字列. 最後のステップを読み込みたい時は`end_step`とする.
- `ix0`, `jx0`, `kx0`: それぞれの方向のMPIスレッドの数
- `ix00`, `jx00`, `kx00` : マージンを除いたアップグレードした格子点数
- `ununiform_flag` : 鉛直方向に非一様格子点間隔を使うか
- `ix_ununi` : 非一様格子点間隔を使うときの一様な領域での格子点の数
- `dx00` : 非一様格子点間隔を使うときの一様な領域での格子点間隔

計算領域に関わるパラメタのみ`upgrade_param_set.F90`で設定している.

```fortran
  upgd%xmax = 0.94d0*rstar; upgd%xmin = 0.71d0*rstar
  upgd%ymax = orgl%ymax; upgd%ymin = orgl%ymin
  upgd%zmax = orgl%zmax; upgd%zmin = orgl%zmin
```