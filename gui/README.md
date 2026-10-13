# VENUS2026 GUI

桌面工具(PyQt6),两个功能页:

1. **输入文件生成** —— 左侧参数表单(体系定义/采样方案/运行控制/束流/浴谱/成员参数/容器键,
   可折叠分组:体系定义/采样方案/输出恒展开,浴谱/束流/成员参数/容器键默认折叠,
   折叠不禁用控件——已填值仍参与校验与写出;校验出错的键所在折叠组自动展开;
   选中带成员键的方案词自动展开成员参数组),
   右侧起始散射构型预览(3D OpenGL 或 2D 投影后备;平衡几何 + 抽样代表值,示意)。
   表单值语义 = "空 = 不写出该键"(与 Fortran 键表可省默认对齐);校验器移植全部
   写盘前可静态判定的 fatal 规则(必填/类型/槽宽/范式窗口/束流搭配/公度/成员正性/
   积分器×阶域…)。
   条件显隐(2026-10-13 裁决):被其余设置判"不需要"的行整行隐藏并清空
   (范式键/定点落点坐标/自适应步长界/浴参数/谱表文件/平动温度/速度展宽/固定方位角/
   成员内二态分支),隐藏键绝不写出;ORDER 等多选一选项为下拉,积分阶选项随积分器变
   (VERLET 1/2·SYMPLECTIC 4/6/8·RADAU 1-15);行标签显示物理量,KEY 名入悬停提示。
   行布局(2026-10-13 裁决):同类设置一行多列自然成组(碰撞能·初始间距、角色对、
   坐标对、三计数键、积分三连、模式+值、成员按词打包每行≤3 键),隐藏某键时同行
   其余键左移、整行键全隐藏时行塌缩;模式类下拉选项显示"值 · 含义"(写出取原值);
   旁注并入悬停提示。表单不提供主势能面/势能面数/动能·势能·总能输出键
   (2026-10-13 裁决移除;键仍在 Fortran 键表,运行走默认)。
   界面文案约定(2026-10-13 裁决):书面措辞——"表面范式/分子范式"作"表面散射/气相散射"、
   "禁写"作"不支持"、字段必填以行标签**加粗星号**标注(条件必填的星号随范式增减),
   不在占位符里写"必填"字样;可见文案书面简洁,括注细节入悬停(2026-10-13 裁决)。
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

## 主题与样式(2026-10-13 裁决: 方案 A 令牌化浅色主题 + 图标 + 密度重排)

唯一样式入口 = `gui/theme.py`:色彩令牌(角色化命名,交互蓝 `#2563EB` 为强调色)、
间距标尺 `SP_S/M/L/XL = 4/8/12/16`、三档排版(正文 9 / 组头 10 / 说明 8 pt);
`apply_theme()` 一次装配 Fusion 底座 + Microsoft YaHei UI 字体 + 调色板 + 全局
QSS。**页面代码零 `setStyleSheet`**(`test_theme.py` 防回漂看守)——局部差异走
objectName(`#SectionHeader` 分组标题色带)或 dynamic property:

| 标记 | 载体 | 效果 |
| --- | --- | --- |
| `class="primary"` | 写出/开始编译按钮 | 实心强调主按钮(其余描边次级) |
| `class="caption"` | 图例/扫描信息/注记/栏标题 | 弱化说明文字(小号灰) |
| `class="card"` | 预览画布外框 | 卡片(白底圆角描边) |
| `class="log"` | 编译日志 | 等宽字体 |
| `state="ok"/"error"` | 校验结果框+摘要行 | 绿/红语义色(文本不变, 仅样式) |
| `objectName=GlScaleLabel` | 3D 预览刻度角标 | 白底描边小角标(GL 视口右上) |

3D 预览视口(2026-10-13 渲染增强):背景/轴色走令牌 `gl_bg`/`gl_axis`
(浅而不近白, 白系原子可辨), 原子球色查表集中在 `preview_view.py` 的
SCHEMES 方案表(CPK 经典浅底适配 / Okabe-Ito 色盲安全系衍生, 工具条
"配色"下拉切换, 默认 CPK), 带尺度 xyz 轴 + "刻度 N Å"角标。

图标为 `gui/assets/icons/` 手绘 SVG(应用/页签/按钮/折叠 chevron),`theme.icon(名)`
按名取用、缺文件返回空图标不崩。语义红线: 标签文本/tooltip/显隐/收集/折叠契约
不受样式影响(既有测试看守);2D 投影画布三处 RGB 为测试锚,不随主题变。

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

## 单位制(2026-10-13 起随 master)

内部单位制为 **eV/Å/amu/fs**(2026-10-06 裁决方案 A,2026-10-13 并入 master):
输入面 E_REL/E_VIB/E_STAB/E_BAR/E_SCALE 按 [eV],DT/DT_BATH/DT_EQ 按 [fs],
速度 [Å/fs](v = √(2·E_REL/(μ·dyn_ev)),dyn_ev = 103.6426965,`gui/core/preview.py`
镜像 `src/globals/consts.f90` 且有镜像测试防漂移)。旧 kcal/mol+(10 fs) 制的
读数归档见 `cases/case_spectrum_source/baseline_kcal_v1/`。

## 防漂移约定

- `keyspec.py` ↔ `src/input.f90` 键表/枚举/分派表:`test_keyspec.py` 逐行比对源文本。
- 元素质量表 ↔ `specio.f90` DATA 表:逐值比对(程序物理常数,不得重推)。
- `pes_scan.py` ↔ `scripts/prebuild.sh`:`test_pes_scan.py::test_prebuild_differential`
  在临时 systems/ 树上同时跑两者,判定(通过/中止)必须一致。

## 注意

- 窗口适配(2026-10-13):初始尺寸钳制在屏幕可用几何内(下沿不压任务栏),两页
  内容各处滚动区(最小高 236px 起步);左栏内容最小宽 502-568px(主题批间距后
  实测,gui/logs/theme_shots.log)< 默认分栏份额约 687px,平时不出横向滚动条
  (用户拉窄分栏时滚动条仍作后备)。
- 预览的抽样量(N_B/N_AIM/N_THTA/N_CHI 采样模式)以代表值代入,视图下方有
  "示意"注记;片段内部几何取体系文件平衡坐标,不含振动/转动采样。
- 预览默认 3D(OpenGL);个别显卡/驱动不兼容时切换"2D 投影"后备。
- 编译验证的隔离沙盒默认保留(失败诊断 + build.log 落盘);"清理沙盒"按钮删除。
