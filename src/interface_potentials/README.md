# interface_potentials/ — 甲类力成员（势能面扩展点）

接口 A 注册族成员：每成员一个 `pes_<名>.f90` 文件，自带 `params_t`
（成员=代码+数据），力例程按甲类签名交付 **eV 力核**（g_add[eV/Å],
v_add[eV]）；×23.0605×C1 换算集中在 facade（单位契约，主计划 §5 N1）。

接入 = 新文件 + 装配相位一行 `call facade_reg_pes('<名>', <例程>)`，
既有代码零触及（T2 判据演练点）。

C2 计划成员（算法自旧库移植，只读参考）：pes_reference（HARMONIC/MORSE/LEPS）、
pes_emt_nn、pes_rst、pes_glo_spring、pes_langevin_ghost、pes_ldfa。
参考质料：interface_TEST.f90 / EMT-NN / RST / VERLET.f90 内嵌力。
