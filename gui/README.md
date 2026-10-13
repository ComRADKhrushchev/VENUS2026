# VENUS2026 GUI

桌面工具(PyQt6),两个功能页:

1. **输入文件生成** —— 左侧参数表单(运行控制/体系定义/采样方案/束流/浴谱/成员参数/容器键),
   右侧起始散射构型预览(3D OpenGL 或 2D 投影后备;平衡几何 + 抽样代表值,示意)。
   表单值语义 = "空 = 不写出该键"(与 Fortran 键表可省默认对齐);校验器移植全部
   写盘前可静态判定的 fatal 规则(必填/类型/槽宽/范式窗口/束流搭配/公度/成员正性…)。
   界面文案约定(2026-10-13 裁决):书面措辞——"表面范式/分子范式"作"表面散射/气相散射"、
   "禁写"作"不支持"、字段必填以行标签**加粗星号**标注(条件必填的星号随范式增减),
   不在占位符里写"必填"字样。
2. **PES 容器检查** —— 对自定义 Tier 1 纯势容器文件夹:①十名扫描 + 八中止
   (prebuild.sh 语义的 Python 镜像);②缺失可选槽位按设定生成 Fortran 代码
   (`gui_slots.f90`:term 三变体 / classify 两变体 / keys+load 接线;方法组五名不生成);
   ③隔离树编译验证(src+scripts+reference 临时副本里激活容器 make,绝不触碰工作树;
   reference 与被检容器并存对编译合法)。

## 启动

```
cd <仓库根>
python -m pip install -r gui/requirements.txt
python -m gui.main          # 或 python gui/main.py
```

编译验证需要 Git Bash + make + gfortran(ifx 可选)。bash 定位优先
`C:\Program Files\Git\...`,PATH 命中的 WSL bash 会被排除;可用环境变量
`VENUS_GUI_BASH` 显式指定。

## 结构

```
gui/main.py            入口(QMainWindow + 两页 QTabWidget)
gui/pages/             input_page / pes_page(Qt 界面层)
gui/view3d/            preview_view(3D GL 主 + QPainter 2D 投影后备)
gui/core/              无 Qt 依赖的逻辑核(可独立单测):
  keyspec.py           键表+枚举+成员分派(镜像 src/input.f90;同步测试防漂移)
  sysdir.py            system/ 解析+范式推导(镜像 utils/specio.f90 + sysdef_interface)
  inputwrite.py        input_qct.txt 生成器 + 语法镜像解析(金标:在档案例)
  validate.py          跨键静态校验(全部写盘前可判 fatal 规则)
  preview.py           初始构型装配(beam_laws/inc_surface/inc_pair 公式移植)
  pes_scan.py          十名扫描+八中止(prebuild.sh:85-168 镜像;差分测试对拍)
  pes_gen.py           缺失槽位代码生成器(term/classify/keys+load 模板)
  sandbox.py           隔离树编译(临时副本 make,日志落盘)
gui/tests/             pytest(常规批默认;`-m slow` 为整树编译慢测)
```

## 测试

```
python -m pytest gui/tests -q        # 常规(含键表/质量表镜像同步、案例金标、装配手算、prebuild 差分、offscreen 冒烟)
python -m pytest gui/tests -m slow -q  # 慢测(auto_sample 隔离 make / 生成槽位端到端编译)
```

## 防漂移约定

- `keyspec.py` ↔ `src/input.f90` 键表/枚举/分派表:`test_keyspec.py` 逐行比对源文本。
- 元素质量表 ↔ `specio.f90` DATA 表:逐值比对(程序物理常数,不得重推)。
- `pes_scan.py` ↔ `scripts/prebuild.sh`:`test_pes_scan.py::test_prebuild_differential`
  在临时 systems/ 树上同时跑两者,判定(通过/中止)必须一致。

## 注意

- 预览的抽样量(N_B/N_AIM/N_THTA/N_CHI 采样模式)以代表值代入,视图下方有
  "示意"注记;片段内部几何取体系文件平衡坐标,不含振动/转动采样。
- 预览默认 3D(OpenGL);个别显卡/驱动不兼容时切换"2D 投影"后备。
- 编译验证的隔离沙盒默认保留(失败诊断 + build.log 落盘);"清理沙盒"按钮删除。
