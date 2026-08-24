# SCUFF-EM NFRHT Fin Structure：Poynting Vector Distribution 后处理与可视化流程设计文档

本文档用于指导 AI coder 实现一个面向 SCUFF-EM `scuff-neq` 结果的后处理脚本/小型 Python 项目。目标是从多频率分片任务目录中读取并拼合 spatially-resolved flux 结果，生成与原始输出目录结构解耦的明文中间数据集，并在此基础上进行频率积分和 Poynting vector 分布可视化。

本文档强调“流程设计、数据契约、错误处理与可视化约束”，不提供可直接照抄的完整代码实现。实现者应根据本文档建立清晰模块边界，避免把 SCUFF 输出读取、频率拼合、物理后处理和绘图逻辑耦合在一起。

---

## 1. 背景与目标

用户使用 SCUFF-EM 模拟 3 周期 fin structure 的 near-field radiative heat transfer（NFRHT）。几何结构本质上是 `xy` 平面中的二维 interpenetrating slabs。为了在 SCUFF-EM 中计算，几何体在 `z` 方向上被人为赋予厚度，例如：

```text
zmin = 0
zmax = 1.5
```

`EPFile` 指定了多个关键截面上的 evaluation points，例如：

```text
z = 0.75  # detailed
z = 1.0   # fewer evpoints
z = 1.5   # fewer evpoints
```

计算结果目录具有如下结构：

```text
task/
  subtask-c1-1/
    evpoints/
    mshs/
    logs/
    omegas/
    run.sbatch
    scuff-neq.log
    scuff-analyze.log
    task.pp
    task.scuffgeo
    task.SRFlux
    task.SIFlux.EMTPFT
    silica-14p-sub1.dat

  subtask-c1-2/
    ...

  subtask-c1-3/
    ...
```

其中：

```text
c1 表示同一组配置；
1, 2, 3, ... 表示该配置下为了提高计算效率而拆分的 frequency-array division index。
```

后处理程序需要完成：

1. 根据必须指定的 `--group-name=c1` 自动发现并读取 `subtask-c1-1/2/3/...`。
2. 读取各分片中的 `task.SRFlux`，并拼合频率分片。
3. 将拼合后的数据转成与原始目录结构无关的明文 canonical dataset。
4. 可选地调用 `scuff-integrate` 或其他后端进行频率积分。
5. 默认输出 integrated Poynting vector over frequency domain 的二维截面分布图。
6. 可视化阶段必须适应 evaluation points 的非完整矩形区域特征，即“大矩形区域中挖去若干小矩形孔洞”的采样网格。

---

## 2. 设计总原则

本项目必须采用分层设计。推荐架构如下：

```text
SCUFF raw output
  ↓
Reader / Collector
  只关心 SCUFF 输出目录结构和文件格式
  ↓
Canonical plain-text dataset
  与 subtask-* 目录结构解耦
  ↓
Postprocessor
  做频率积分、温度加权、source/transform 汇总等物理后处理
  ↓
Visualization
  只读取后处理后的 canonical dataset，不直接读取 SCUFF 原始输出
```

### 2.1 关键约束

实现时必须遵守以下约束：

```text
1. 读取层和后处理层彻底解耦。
2. 后处理层不得依赖 subtask-c1-1/2/3 目录结构。
3. 可视化层不得直接读取 task.SRFlux。
4. 中间数据禁止使用二进制格式。
5. 中间数据必须是可检查、可 diff、可长期保存的明文格式。
6. 频率分片拼合不得简单拼接文本文件。
7. 可视化不得默认使用 contourf、tricontourf 或 griddata 插值。
8. 缺失 evaluation points 对应区域必须保持空白或 mask，不得插值补洞。
```

推荐使用：

```text
TSV / CSV / YAML / JSONL
```

禁止使用：

```text
pickle / npy / npz / HDF5 / parquet / feather / 自定义 bin
```

---

## 3. 推荐项目结构

建议项目结构如下：

```text
scuff_pv_post/
  pyproject.toml
  README.md

  scuff_pv_post/
    __init__.py

    cli.py
      # 命令行入口

    discovery.py
      # 发现 subtask-group-index 目录

    srflux_reader.py
      # 读取和解析 *.SRFlux

    column_mapping.py
      # 解析 header / column_map.yaml / fallback 列映射

    canonical_dataset.py
      # 生成 merged_srflux.tsv, points.tsv, omegas.tsv, manifest.yaml

    validation.py
      # 完整性检查、重复检查、坐标网格检查、频率检查

    integration.py
      # 调用 scuff-integrate 或 Python fallback

    pvmst_reader.py
      # 读取 integrated PVMST 或自定义 integrated TSV

    visualization.py
      # masked cell map / quiver / scatter fallback

    report.py
      # 生成 error report / data summary / provenance

  examples/
    column_map_srflux.yaml
    plot_config.yaml
    TemperatureFile.example

  tests/
    test_discovery.py
    test_merge_keys.py
    test_masked_grid.py
    test_duplicate_detection.py
```

不要求完全照此命名，但模块职责必须等价清晰。

---

## 4. 命令行接口设计

建议提供三个主要子命令：

```bash
scuff-pv collect
scuff-pv integrate
scuff-pv plot
```

也可以提供一个总流程命令：

```bash
scuff-pv run-all
```

但 `run-all` 只能作为 wrapper，内部仍然调用分层模块，不得把全部逻辑写在一个函数里。

---

## 5. Collect 阶段：读取与拼合频率分片

### 5.1 输入

典型命令：

```bash
scuff-pv collect \
  --task-dir ./task \
  --group-name c1 \
  --filebase task \
  --out ./pv_out/c1
```

参数含义：

| 参数 | 必需 | 含义 |
|---|---:|---|
| `--task-dir` | 是 | 包含 `subtask-*` 子目录的根目录 |
| `--group-name` | 是 | 配置组名称，例如 `c1` |
| `--filebase` | 否 | 默认 `task`，用于寻找 `task.SRFlux` |
| `--column-map` | 否 | 显式列映射 YAML |
| `--out` | 是 | canonical dataset 输出目录 |
| `--coord-tol` | 否 | 坐标归一化容差，用于生成稳定 point_id |
| `--duplicate-policy` | 否 | 重复 key 处理策略，默认 `error` |

### 5.2 子目录发现规则

给定：

```text
--group-name c1
```

只允许匹配：

```text
subtask-c1-1
subtask-c1-2
subtask-c1-3
...
```

不得匹配：

```text
subtask-c10-1
subtask-c1-extra
subtask-c2-1
subtask-c1
```

推荐使用严格正则：

```text
^subtask-{group_name}-(\d+)$
```

其中最后的数字为 frequency division index，必须用于排序。

### 5.3 每个 subtask 中的文件选择规则

默认优先读取：

```text
{subtask_dir}/{filebase}.SRFlux
```

例如：

```text
subtask-c1-1/task.SRFlux
```

如果该文件不存在，可以 fallback 到：

```text
*.SRFlux
```

但只有在唯一匹配时才允许自动选择。

如果存在多个 `*.SRFlux` 且无法由 `--filebase` 唯一确定，必须报错，不得随机选择。

### 5.4 不建议读取的文件

`task.SIFlux.EMTPFT` 是 spatially-integrated flux，可作为全局热流校验，但不是局域 Poynting vector 分布的主要输入。

`mshs/` 和 `task.scuffgeo` 可用于未来叠加几何边界，但不应作为最小流程的强依赖。

`scuff-neq.log` 和 `scuff-analyze.log` 可用于 provenance、debug 和错误报告，但不应影响正常数据读取流程，除非用户显式启用 log consistency check。

---

## 6. SRFlux 解析策略

### 6.1 不允许硬编码列号

不同 SCUFF-EM 版本、不同输出选项或不同 flux 类型可能导致列数或列含义发生变化。实现时不得假设“第 9 列永远是 Sz”。

推荐解析优先级：

```text
1. 从 SRFlux 文件 header 自动识别列名。
2. 如果 header 不足，读取用户提供的 column_map.yaml。
3. 如果二者都不可用，使用内置 fallback 映射，但必须显示 warning，并在 manifest.yaml 中记录。
```

### 6.2 column_map.yaml 的设计

示例：

```yaml
format: scuff_neq_srflux
index_base: 1

columns:
  omega: 1
  transform_id: 2
  source_id: 3
  x: 4
  y: 5
  z: 6
  Sx_flux: 7
  Sy_flux: 8
  Sz_flux: 9
  Mxx_flux: 10
  Mxy_flux: 11
  Mxz_flux: 12
  Myx_flux: 13
  Myy_flux: 14
  Myz_flux: 15
  Mzx_flux: 16
  Mzy_flux: 17
  Mzz_flux: 18
```

说明：

- `index_base=1` 表示列号按人类习惯从 1 开始。
- 实现内部可以转为 0-based index。
- 如果某些 Maxwell stress tensor 列不存在，不应影响 Poynting vector 可视化；但必须在 manifest 中记录缺失列。

### 6.3 header 保留

读取 `.SRFlux` 时，应保留 header 信息到：

```text
manifest.yaml
raw_headers/
```

推荐输出：

```text
pv_out/c1/raw_headers/subtask-c1-1.SRFlux.header.txt
pv_out/c1/raw_headers/subtask-c1-2.SRFlux.header.txt
```

目的：

1. 便于回溯 SCUFF 输出格式。
2. 便于人工检查列映射是否正确。
3. 便于 debugging。

---

## 7. Canonical dataset 设计

Collect 阶段必须输出一个独立的 canonical dataset。后处理和可视化只能依赖该 dataset。

推荐输出：

```text
pv_out/c1/
  manifest.yaml
  merged_srflux.tsv
  points.tsv
  omegas.tsv
  sources.tsv
  transforms.tsv
  duplicate_report.tsv
  missing_report.tsv
  raw_headers/
```

其中只有前三个是最小必需：

```text
manifest.yaml
merged_srflux.tsv
points.tsv
omegas.tsv
```

---

## 8. 点坐标索引设计

### 8.1 不要依赖 EPFile 行号

同一组 computation 的不同 frequency subtask 中，evaluation points 理论上应一致，但不能假设行号永远一致。

必须使用坐标生成稳定 `point_id`。

推荐逻辑：

```text
point_id = hash(round(x / coord_tol), round(y / coord_tol), round(z / coord_tol))
```

`coord_tol` 默认可设为：

```text
1e-9 或 1e-8
```

实际取值应允许用户通过命令行配置。

### 8.2 points.tsv

推荐格式：

```text
point_id    x          y          z          plane_id     section_label
p_xxx       -1.0000    -0.5000    0.7500     z_0.750      mid_z
p_yyy       -0.9900    -0.5000    0.7500     z_0.750      mid_z
p_zzz       -1.0000    -0.5000    1.0000     z_1.000      upper_inner
```

其中：

- `point_id` 是稳定索引。
- `plane_id` 根据 `z` 坐标生成。
- `section_label` 可选，可由用户配置。

### 8.3 坐标容差与重复点

如果不同原始文件中出现坐标极其接近但不完全相同的点，应按 `coord_tol` 归一化后判断是否为同一点。

如果多个点归一化到同一个 `point_id`，但原始坐标偏差超过 `coord_tol` 的合理范围，应报 warning 或 error。

推荐规则：

```text
同一 point_id 下：
  max(abs(x - x_mean)) <= coord_tol
  max(abs(y - y_mean)) <= coord_tol
  max(abs(z - z_mean)) <= coord_tol
```

否则视为 coordinate collision。

---

## 9. 频率索引设计

### 9.1 omegas.tsv

推荐格式：

```text
omega_id    omega          subtask_index    source_file
w000001     0.010000       1                subtask-c1-1/task.SRFlux
w000002     0.012000       1                subtask-c1-1/task.SRFlux
w000101     0.200000       2                subtask-c1-2/task.SRFlux
```

### 9.2 频率排序

最终 `merged_srflux.tsv` 必须按数值频率排序，而不是按文件顺序排序。

排序 key 推荐：

```text
omega, transform_id, source_id, point_id
```

### 9.3 频率完整性检查

Collect 阶段应检查：

1. 是否存在重复 `omega`。
2. 是否存在非单调频率分片。
3. 是否存在明显频率缺口。
4. 每个 `omega` 下的 point/source/transform 组合是否完整。

注意：频率分片之间允许边界频率重复，但重复数据必须一致，或者按明确策略处理。

---

## 10. merged_srflux.tsv 设计

推荐核心列：

```text
omega
transform_id
source_id
point_id
x
y
z
Sx_flux
Sy_flux
Sz_flux
Mxx_flux
Mxy_flux
Mxz_flux
Myx_flux
Myy_flux
Myz_flux
Mzx_flux
Mzy_flux
Mzz_flux
subtask_index
source_file
```

最小必需列：

```text
omega
transform_id
source_id
point_id
x
y
z
Sx_flux
Sy_flux
Sz_flux
```

推荐 primary key：

```text
(omega, transform_id, source_id, point_id)
```

若后续需要区分更多维度，例如 polarization、body pair、temperature configuration，可扩展 key，但不得破坏已有字段含义。

---

## 11. 重复数据处理

### 11.1 重复 key 的定义

重复 key 指：

```text
(omega, transform_id, source_id, point_id)
```

完全相同的多行。

### 11.2 默认策略

默认策略必须是：

```text
--duplicate-policy error
```

即：

- 如果重复行所有数值列完全一致：可保留一条，同时记录到 `duplicate_report.tsv`。
- 如果重复行数值列不一致：必须报错，不得悄悄平均。
- 只有用户显式指定 `--duplicate-policy mean` 时才允许平均。
- 如果指定 `mean`，必须在 manifest 中记录，并输出 averaged duplicate report。

### 11.3 duplicate_report.tsv

至少包含：

```text
omega
transform_id
source_id
point_id
n_duplicates
source_files
max_abs_diff
action
```

`action` 可取：

```text
kept_first
averaged
error
ignored_identical
```

---

## 12. 缺失数据检查

### 12.1 允许的缺失

evaluation points 本身不是完整矩形，而是“大矩形减去若干小矩形孔洞”。这些孔洞在可视化中应保持 mask，不应视为错误。

### 12.2 不允许的缺失

对于同一个 `omega`，如果某个 `point_id` 在其他频率都存在，但该频率缺失，应视为潜在错误。

例如：

```text
omega = 0.10 有 10000 个点
omega = 0.12 只有 7300 个点
```

除非用户显式指定该频率采用了不同 EPFile，否则应报 error 或 warning。

### 12.3 missing_report.tsv

建议输出：

```text
omega
transform_id
source_id
expected_n_points
actual_n_points
missing_point_count
action
```

如果存在缺失 point，额外输出：

```text
missing_points_by_omega.tsv
```

包含：

```text
omega
transform_id
source_id
point_id
x
y
z
```

---

## 13. manifest.yaml 设计

`manifest.yaml` 是数据集的核心说明文件，应记录：

```yaml
dataset_version: 1

task_dir: /absolute/or/relative/path/to/task
group_name: c1
filebase: task

created_at: "YYYY-MM-DDTHH:MM:SS"
software:
  project_name: scuff_pv_post
  project_version: "0.1.0"
  python_version: "..."
  scuff_integrate_path: "..."
  scuff_integrate_version: "unknown or parsed"

discovery:
  matched_subtasks:
    - subtask-c1-1
    - subtask-c1-2
  ignored_subtasks:
    - subtask-c2-1

column_mapping:
  method: header | user_yaml | fallback
  column_map_file: column_map_srflux.yaml
  warning: null

coordinate:
  coord_tol: 1e-9
  point_id_method: rounded_coordinate_hash

merge:
  primary_key:
    - omega
    - transform_id
    - source_id
    - point_id
  duplicate_policy: error
  n_rows_raw: 1000000
  n_rows_merged: 999990
  n_points: 10000
  n_omegas: 100
  n_sources: 2
  n_transforms: 1

validation:
  status: pass | warning | fail
  reports:
    duplicate_report: duplicate_report.tsv
    missing_report: missing_report.tsv

integration:
  status: not_run | success | failed
  backend: scuff-integrate | python
  temperature_file: TemperatureFile
  output_file: merged.PVMST

visualization:
  default_mode: masked_cellmap_quiver
  interpolation: forbidden_by_default
```

---

## 14. Integration 阶段：频率积分

### 14.1 物理含义

`.SRFlux` 是 frequency-resolved spatially-resolved flux。若目标是 integrated Poynting vector over frequency domain，需要进行频率积分。对于 NFRHT 问题，频率积分通常还涉及温度配置和热占据因子。

因此，默认不应把 `.SRFlux` 的数值直接对频率求和并称为最终热流。推荐默认使用 SCUFF-EM 提供的 `scuff-integrate`。

### 14.2 推荐命令

```bash
scuff-pv integrate \
  --dataset ./pv_out/c1 \
  --temperature-file ./TemperatureFile \
  --backend scuff-integrate
```

### 14.3 输入与输出

输入：

```text
pv_out/c1/merged.SRFlux
TemperatureFile
```

输出：

```text
pv_out/c1/merged.PVMST
pv_out/c1/pvmst_integrated.tsv
pv_out/c1/integration_report.yaml
```

### 14.4 关于 merged.SRFlux

虽然 canonical dataset 的主文件是：

```text
merged_srflux.tsv
```

但为了调用 `scuff-integrate`，需要额外导出一个 SCUFF-compatible 的：

```text
merged.SRFlux
```

该文件应尽量保持 SCUFF 原始 `.SRFlux` 的列顺序和 header 风格，使 `scuff-integrate` 能直接识别。

这意味着 collect 阶段最好同时保存：

1. 语义化的 `merged_srflux.tsv`。
2. SCUFF-compatible 的 `merged.SRFlux`。

两者内容应等价，但用途不同：

| 文件 | 用途 |
|---|---|
| `merged_srflux.tsv` | 稳定、可读、供 Python 后处理 |
| `merged.SRFlux` | 供 `scuff-integrate` 调用 |

### 14.5 scuff-integrate 错误处理

如果 `scuff-integrate` 不存在：

```text
error: SCUFF_INTEGRATE_NOT_FOUND
suggestion: 确认 scuff-integrate 在 PATH 中，或通过 --scuff-integrate 指定路径。
```

如果 `TemperatureFile` 不存在：

```text
error: TEMPERATURE_FILE_NOT_FOUND
```

如果 `scuff-integrate` 返回非零退出码：

```text
error: SCUFF_INTEGRATE_FAILED
保存 stdout/stderr 到 integration_report.yaml 或 logs/scuff-integrate.stderr.txt
```

如果 `merged.PVMST` 没有生成：

```text
error: PVMST_OUTPUT_MISSING
```

---

## 15. PVMST 读取与 integrated dataset

### 15.1 pvmst_integrated.tsv

积分后应读取 `.PVMST` 并转换为语义化明文表：

```text
pvmst_integrated.tsv
```

推荐列：

```text
transform_id
source_id
point_id
x
y
z
Sx
Sy
Sz
Mxx
Mxy
Mxz
Myx
Myy
Myz
Mzx
Mzy
Mzz
```

最小必需列：

```text
transform_id
source_id
point_id
x
y
z
Sx
Sy
Sz
```

### 15.2 source 汇总策略

Poynting vector 的最终可视化是否需要对 source object 求和，应由配置控制，不能硬编码。

推荐参数：

```text
--source-mode separate | sum | selected
```

含义：

| 模式 | 含义 |
|---|---|
| `separate` | 每个 source 单独输出图 |
| `sum` | 对所有 source 求和后输出总 Poynting vector |
| `selected` | 只选指定 source_id |

默认建议：

```text
--source-mode sum
```

但必须在图标题和文件名中明确标明：

```text
source=sum
```

### 15.3 transform 汇总策略

同理，transform 也不应硬编码。

推荐参数：

```text
--transform-mode separate | sum | selected
```

如果通常只有一个 transform，可默认 `selected` 或 `separate`，但仍需保留扩展性。

---

## 16. Visualization 阶段：适应挖孔矩形采样区域

### 16.1 数据特点

用户的 evaluation points 不是完整矩形区域，而是：

```text
大矩形区域中挖去若干小矩形区域后得到的剩余区域均匀取点。
```

这意味着：

- `x, y` 整体看起来像规则网格。
- 但某些矩形区域内部没有采样点。
- 如果直接使用 `contourf`、`tricontourf`、`griddata` 或其他插值方法，会跨越孔洞产生虚假的连续场。
- 这会导致图像看起来更平滑，但物理上不可信。

因此默认可视化必须使用 sampled-field preserving 策略。

---

## 17. 默认绘图策略：masked cell map

### 17.1 基本原则

默认 colored figure 必须使用：

```text
masked cell map
```

即：

```text
1. 根据真实 x_unique, y_unique 重建规则网格。
2. 真实存在 evpoint 的位置填入数值。
3. 不存在 evpoint 的位置填入 NaN。
4. 用 pcolormesh 或 imshow(nearest) 绘制。
5. NaN 区域保持空白或透明。
```

禁止默认使用：

```text
contourf
tricontourf
griddata
Delaunay triangulation interpolation
streamplot
```

其中 `streamplot` 也不推荐默认使用，因为它默认假设矢量场定义在连续矩形区域上，对挖孔区域容易产生误导。

### 17.2 默认图类型

每个 z-plane 默认输出：

```text
1. cellmap_absSxy.png
   color = |Sxy| = sqrt(Sx^2 + Sy^2)

2. cellmap_Sz.png
   color = Sz

3. cellmap_absSxy_quiver.png
   color = |Sxy|，叠加 Sx,Sy 箭头

4. cellmap_log_absSxy_quiver.png
   color = log10(|Sxy| + eps)，叠加 Sx,Sy 箭头
```

推荐主图：

```text
cellmap_absSxy_quiver.png
```

因为它同时呈现：

```text
color: 局域横向 Poynting vector 强度
arrow: xy 平面内能流方向
```

### 17.3 默认截面

对该类人为挤出厚度的二维结构，默认主截面应为：

```text
z = 0.75
```

因为它是 `zmin=0, zmax=1.5` 的中截面，最能代表二维结构内部的能流分布。

辅助截面：

```text
z = 1.0
z = 1.5
```

其中 `z=1.5` 位于上边界，应谨慎解释，可能包含边界截断效应或表面附近效应。

---

## 18. Grid reconstruction 细节

### 18.1 判断规则网格

对单个 z-plane：

```text
x_unique = sorted(unique(x))
y_unique = sorted(unique(y))
```

如果点集来自规则网格挖孔，则大多数点应满足：

```text
x ∈ x_unique
y ∈ y_unique
```

并且相邻 x/y 间距大致稳定。

应检查：

```text
dx_values = diff(x_unique)
dy_values = diff(y_unique)
```

如果：

```text
max(abs(dx_values - median(dx_values))) <= grid_tol
max(abs(dy_values - median(dy_values))) <= grid_tol
```

则可使用 `pcolormesh`。

如果不是规则网格，则 fallback 到 scatter cell map。

### 18.2 构建网格

对于规则网格：

```text
value_grid[ny, nx] = NaN
Sx_grid[ny, nx] = NaN
Sy_grid[ny, nx] = NaN
Sz_grid[ny, nx] = NaN
```

对每个 point：

```text
ix = index of x in x_unique
iy = index of y in y_unique

value_grid[iy, ix] = value
Sx_grid[iy, ix] = Sx
Sy_grid[iy, ix] = Sy
Sz_grid[iy, ix] = Sz
```

没有 point 的位置保持 NaN。

### 18.3 pcolormesh cell edges

对于等间距网格：

```text
dx = median(diff(x_unique))
dy = median(diff(y_unique))

x_edges = [x0 - dx/2, x1 - dx/2, ..., x_last + dx/2]
y_edges = [y0 - dy/2, y1 - dy/2, ..., y_last + dy/2]
```

如果间距不完全等间距，可根据相邻点中点构造 edges，但若偏离太大应 fallback 到 scatter。

---

## 19. Quiver 叠加策略

### 19.1 箭头向量

默认箭头：

```text
U = Sx
V = Sy
```

箭头只画在：

```text
isfinite(value) and isfinite(Sx) and isfinite(Sy)
```

的位置。

禁止在 NaN / mask 区域画箭头。

### 19.2 箭头稀疏化

由于 evaluation points 可能非常密集，必须支持：

```text
--stride 3
```

或：

```text
--max-arrows 2000
```

推荐默认：

```text
stride = 3
```

但如果点数很少，应自动降低 stride。

### 19.3 箭头归一化

应提供两种模式：

```text
--arrow-mode physical
--arrow-mode normalized
```

含义：

| 模式 | 含义 |
|---|---|
| `physical` | 箭头长度与 Poynting vector 大小成比例 |
| `normalized` | 箭头只表示方向，长度统一或弱依赖大小 |

默认建议：

```text
physical
```

但对于动态范围极大的数据，`normalized` 更适合看流向。

### 19.4 color 与 arrow 的分工

推荐：

```text
color 表示 magnitude 或指定标量；
arrow 表示方向。
```

不要把箭头颜色也映射成另一个复杂变量，避免图像难以解释。

---

## 20. Color value 设计

Visualization 阶段应支持以下 scalar：

| 名称 | 定义 |
|---|---|
| `Sx` | x 方向 Poynting vector |
| `Sy` | y 方向 Poynting vector |
| `Sz` | z 方向 Poynting vector |
| `absSxy` | `sqrt(Sx^2 + Sy^2)` |
| `absS` | `sqrt(Sx^2 + Sy^2 + Sz^2)` |
| `log_absSxy` | `log10(absSxy + eps)` |
| `log_absS` | `log10(absS + eps)` |

默认：

```text
--color-by absSxy
```

### 20.1 signed quantity

对于 `Sx`, `Sy`, `Sz` 这类 signed quantity：

```text
colorbar 应使用 diverging colormap，且推荐以 0 为中心。
```

但具体 colormap 不应强制写死，可以配置。

### 20.2 magnitude quantity

对于 `absSxy`, `absS`, `log_absSxy`：

```text
colorbar 使用 sequential colormap。
```

### 20.3 robust clipping

强烈建议支持：

```text
--vmin-percentile 1
--vmax-percentile 99
```

用于避免少数极端点破坏整体色阶。

但必须在图标题或 metadata 中记录 clipping 信息，例如：

```text
color scale clipped to 1st-99th percentile
```

---

## 21. 非规则网格 fallback：scatter cell map

如果某个截面无法可靠重建为规则网格，禁止自动切换到 contour。

fallback 应为：

```text
scatter cell map
```

即：

```text
plt.scatter(x, y, c=value, marker="s")
```

点大小根据估计的局部 spacing 设置。

原则仍然是：

```text
只在真实 evpoint 位置显示 color，不插值，不补洞。
```

---

## 22. 几何边界叠加

第一版可不实现几何边界叠加，但应预留接口。

建议参数：

```text
--geometry-overlay none | scuffgeo | polygon-file
```

第一版默认：

```text
none
```

后续可从 `task.scuffgeo` 或用户提供的 polygon 文件中读取 slab 边界，用黑色线框或灰色 patch 叠加。

注意：几何边界叠加只是辅助解释，不得影响数值图本身。

---

## 23. Plot 命令设计

推荐命令：

```bash
scuff-pv plot \
  --dataset ./pv_out/c1 \
  --planes 0.75,1.0,1.5 \
  --mode cellmap-quiver \
  --color-by absSxy \
  --vector Sx,Sy \
  --source-mode sum \
  --transform-mode separate \
  --no-interpolation \
  --mask-missing \
  --stride 3 \
  --out ./pv_out/c1/figures
```

默认参数应等价于：

```text
--mode cellmap-quiver
--color-by absSxy
--vector Sx,Sy
--source-mode sum
--no-interpolation true
--mask-missing true
--contour false
--streamplot false
```

如果用户显式请求 contour，应要求：

```text
--allow-interpolation
```

并输出 warning：

```text
Contour/interpolation may create artificial fields across holes in the evaluation-point domain.
```

---

## 24. 输出文件命名规则

建议文件名包含：

```text
group
z-plane
source-mode
transform
color quantity
plot mode
```

例如：

```text
figures/
  c1_z0.750_source-sum_transform-0_absSxy_cellmap.png
  c1_z0.750_source-sum_transform-0_absSxy_cellmap_quiver.png
  c1_z0.750_source-sum_transform-0_Sz_cellmap.png
  c1_z0.750_source-sum_transform-0_log_absSxy_cellmap_quiver.png
```

同时输出：

```text
figures/plot_manifest.yaml
```

记录每张图的输入数据、参数、色阶范围、stride、mask 点数等。

---

## 25. Error handling 总体设计

### 25.1 错误等级

推荐错误等级：

| 等级 | 含义 |
|---|---|
| `fatal` | 无法继续执行 |
| `error` | 当前阶段失败 |
| `warning` | 可以继续，但结果需要谨慎 |
| `info` | 普通记录 |

### 25.2 错误报告

每个阶段应输出独立 report：

```text
collect_report.yaml
integration_report.yaml
plot_report.yaml
```

总 report 写入：

```text
manifest.yaml
```

### 25.3 不要静默失败

以下情况不得静默跳过：

1. 没有匹配到任何 subtask。
2. 某个 subtask 缺失 `.SRFlux`。
3. `.SRFlux` 列映射失败。
4. 数值列包含无法解析的非数字。
5. 重复 primary key 但数值不一致。
6. 某些频率点的 evaluation point 数量异常缺失。
7. `scuff-integrate` 失败。
8. plot 阶段指定的 z-plane 没有任何点。
9. masked grid 中所有值均为 NaN。
10. 用户请求 contour 但没有显式允许 interpolation。

---

## 26. Collect 阶段错误处理清单

### 26.1 group-name 相关

情况：

```text
--group-name 未指定
```

处理：

```text
fatal: GROUP_NAME_REQUIRED
```

情况：

```text
未发现 subtask-c1-*
```

处理：

```text
fatal: NO_MATCHING_SUBTASKS
```

情况：

```text
发现 subtask-c1-1, subtask-c1-3，但缺少 subtask-c1-2
```

处理：

```text
warning 或 error，由 --require-contiguous-subtasks 控制。
```

默认建议 warning，因为用户可能故意删除失败任务，但 report 中必须记录。

### 26.2 文件相关

情况：

```text
task.SRFlux 不存在，但目录中存在唯一 *.SRFlux
```

处理：

```text
warning: SRFLUX_FILEBASE_FALLBACK
```

情况：

```text
存在多个 *.SRFlux，无法判断
```

处理：

```text
fatal: AMBIGUOUS_SRFLUX_FILES
```

情况：

```text
SRFlux 文件为空
```

处理：

```text
fatal: EMPTY_SRFLUX
```

### 26.3 列映射相关

情况：

```text
header 无法解析，且没有 column_map.yaml
```

处理：

```text
warning: USING_FALLBACK_COLUMN_MAP
```

情况：

```text
fallback 映射列数超过实际列数
```

处理：

```text
fatal: COLUMN_MAP_OUT_OF_RANGE
```

情况：

```text
缺失 Sx/Sy/Sz
```

处理：

```text
fatal: POYNTING_COLUMNS_MISSING
```

如果只缺 Maxwell stress tensor 列：

```text
warning: MST_COLUMNS_MISSING
```

### 26.4 数据完整性相关

情况：

```text
某些行无法解析为数字
```

处理：

```text
fatal: NUMERIC_PARSE_FAILED
```

情况：

```text
某些关键列存在 NaN
```

关键列包括：

```text
omega, transform_id, source_id, x, y, z, Sx_flux, Sy_flux, Sz_flux
```

处理：

```text
fatal: NAN_IN_REQUIRED_COLUMNS
```

情况：

```text
同一 point_id 出现 coordinate collision
```

处理：

```text
fatal 或 warning，由 --coordinate-collision-policy 控制。默认 fatal。
```

---

## 27. Integration 阶段错误处理清单

情况：

```text
merged.SRFlux 不存在
```

处理：

```text
fatal: MERGED_SRFLUX_NOT_FOUND
```

情况：

```text
TemperatureFile 不存在
```

处理：

```text
fatal: TEMPERATURE_FILE_NOT_FOUND
```

情况：

```text
scuff-integrate 不在 PATH
```

处理：

```text
fatal: SCUFF_INTEGRATE_NOT_FOUND
```

情况：

```text
scuff-integrate 返回非零 code
```

处理：

```text
fatal: SCUFF_INTEGRATE_FAILED
```

并保存：

```text
stdout
stderr
command
return_code
```

情况：

```text
PVMST 文件生成但为空
```

处理：

```text
fatal: EMPTY_PVMST
```

情况：

```text
PVMST 行数与 points/source/transform 组合明显不一致
```

处理：

```text
warning 或 error；默认 error。
```

---

## 28. Plot 阶段错误处理清单

情况：

```text
pvmst_integrated.tsv 不存在
```

处理：

```text
fatal: INTEGRATED_DATASET_NOT_FOUND
```

情况：

```text
指定 z=0.75，但找不到该截面
```

处理：

```text
fatal: PLANE_NOT_FOUND
```

可提供建议：

```text
available planes: 0.750000, 1.000000, 1.500000
```

情况：

```text
指定 color-by 字段不存在
```

处理：

```text
fatal: UNKNOWN_COLOR_QUANTITY
```

情况：

```text
该 plane 中所有 color value 为 NaN
```

处理：

```text
fatal: ALL_VALUES_NAN
```

情况：

```text
网格不是规则网格
```

处理：

```text
warning: IRREGULAR_GRID_FALLBACK_TO_SCATTER
```

不得 fallback 到 contour。

情况：

```text
用户请求 contour，但未指定 --allow-interpolation
```

处理：

```text
fatal: INTERPOLATION_NOT_ALLOWED
```

情况：

```text
quiver arrow 太多
```

处理：

```text
warning: QUIVER_DOWNSAMPLED
```

自动增大 stride 或使用 `--max-arrows` 控制。

---

## 29. 测试要求

实现者至少应写以下测试。

### 29.1 discovery 测试

输入目录：

```text
subtask-c1-1
subtask-c1-2
subtask-c2-1
subtask-c10-1
```

给定：

```text
group_name = c1
```

只应返回：

```text
subtask-c1-1
subtask-c1-2
```

### 29.2 duplicate key 测试

构造两条相同 key：

```text
(omega, transform_id, source_id, point_id)
```

测试：

1. 数值一致时保留一条并记录 duplicate。
2. 数值不一致时默认报错。
3. `--duplicate-policy mean` 时输出平均值并记录。

### 29.3 masked grid 测试

构造：

```text
5 x 5 网格
挖掉中间 2 x 2 区域
```

要求：

1. 输出 grid shape 为 `5 x 5`。
2. 被挖掉区域为 NaN。
3. plot 函数不得填补 NaN。
4. quiver 不得在 NaN 区域画箭头。

### 29.4 irregular grid fallback 测试

构造非规则 x/y 点集，要求：

```text
plot mode 自动 fallback 到 scatter cell map
```

不得调用 contour / interpolation。

### 29.5 z-plane 容差测试

如果 z 坐标存在微小浮点误差：

```text
0.7500000001
0.7499999999
```

应能按 `plane_tol` 归为同一 plane。

### 29.6 integration failure 测试

模拟 `scuff-integrate` 返回非零状态码，要求：

1. 报出 `SCUFF_INTEGRATE_FAILED`。
2. 保存 stderr。
3. 不生成假结果文件。

---

## 30. 推荐实现顺序

建议 AI coder 按以下顺序实现：

```text
1. 实现 discovery：根据 --group-name 找 subtask。
2. 实现 SRFlux 读取：先用 column_map.yaml，不急于做 header 自动解析。
3. 实现 point_id 和 canonical dataset 输出。
4. 实现 duplicate / missing / coordinate validation。
5. 实现 merged.SRFlux 导出，确保 scuff-integrate 可读。
6. 实现 scuff-integrate wrapper。
7. 实现 PVMST → pvmst_integrated.tsv。
8. 实现 masked cell map。
9. 实现 quiver overlay。
10. 实现 scatter fallback。
11. 最后再做 geometry overlay、header 自动识别和高级图形选项。
```

---

## 31. 最小可用版本的验收标准

第一版不需要支持所有高级功能，但必须满足：

```text
1. 能通过 --group-name=c1 自动读取 subtask-c1-*。
2. 能把多个 task.SRFlux 拼合成 merged_srflux.tsv。
3. 能生成 points.tsv, omegas.tsv, manifest.yaml。
4. 能检测重复 key 和明显缺失数据。
5. 能导出 scuff-integrate 可读的 merged.SRFlux。
6. 能调用 scuff-integrate 生成 integrated Poynting vector。
7. 能读取 integrated result 并生成 pvmst_integrated.tsv。
8. 能对 z=0.75, 1.0, 1.5 输出 masked cell map。
9. 能输出 colored figure，不使用 contour interpolation。
10. 能叠加 Sx,Sy quiver，且不在孔洞区域画箭头。
11. 所有中间结果都是明文。
```

---

## 32. 关键物理解释要求

脚本输出的 README 或 report 中应提醒用户：

1. `.SRFlux` 是 frequency-resolved spatially-resolved flux，不应直接等同于最终热平均 Poynting vector。
2. integrated Poynting vector 需要频率积分，通常还需要温度配置和热权重。
3. 对二维 `xy` interpenetrating slabs 人为挤出厚度的模型，`z=0.75` 更适合作为主截面。
4. `z=1.5` 位于边界，解释时应谨慎。
5. 颜色图中的空白区域代表没有 evaluation point，不代表 flux 为零。
6. 若使用 log scale，低值区域会被视觉放大，应在图注中说明。
7. 若使用 percentile clipping，图中颜色不是全局绝对范围，应在图注中说明。
8. quiver 箭头的长度是否代表真实大小取决于 `arrow-mode`，必须在图注中说明。

---

## 33. 不应实现的默认行为

以下行为不应作为默认行为：

```text
1. 自动 contourf 平滑。
2. 自动 griddata 插值。
3. 自动 Delaunay triangulation。
4. 自动补齐孔洞。
5. 自动把 NaN 当成 0。
6. 自动把所有 source_id 混在一起且不记录。
7. 自动忽略重复频率。
8. 自动忽略缺失 point。
9. 自动读取错误 group 的 subtask。
10. 自动选择多个 SRFlux 文件中的一个。
```

如果用户显式要求这些行为，必须通过明确参数开启，并在 report 中记录。

---

## 34. 最终流程示例

推荐完整流程：

```bash
scuff-pv collect \
  --task-dir ./task \
  --group-name c1 \
  --filebase task \
  --column-map ./column_map_srflux.yaml \
  --out ./pv_out/c1

scuff-pv integrate \
  --dataset ./pv_out/c1 \
  --temperature-file ./TemperatureFile \
  --backend scuff-integrate

scuff-pv plot \
  --dataset ./pv_out/c1 \
  --planes 0.75,1.0,1.5 \
  --mode cellmap-quiver \
  --color-by absSxy \
  --vector Sx,Sy \
  --source-mode sum \
  --no-interpolation \
  --mask-missing \
  --stride 3 \
  --out ./pv_out/c1/figures
```

输出目录：

```text
pv_out/c1/
  manifest.yaml
  collect_report.yaml
  integration_report.yaml
  plot_report.yaml

  merged_srflux.tsv
  merged.SRFlux
  points.tsv
  omegas.tsv
  sources.tsv
  transforms.tsv

  pvmst_integrated.tsv
  merged.PVMST

  duplicate_report.tsv
  missing_report.tsv

  raw_headers/
    subtask-c1-1.SRFlux.header.txt
    subtask-c1-2.SRFlux.header.txt

  figures/
    plot_manifest.yaml
    c1_z0.750_source-sum_transform-0_absSxy_cellmap.png
    c1_z0.750_source-sum_transform-0_absSxy_cellmap_quiver.png
    c1_z0.750_source-sum_transform-0_Sz_cellmap.png
    c1_z0.750_source-sum_transform-0_log_absSxy_cellmap_quiver.png
    c1_z1.000_source-sum_transform-0_absSxy_cellmap.png
    c1_z1.500_source-sum_transform-0_absSxy_cellmap.png
```

---

## 35. 一句话总结

该后处理项目的核心不是“把 SCUFF 输出读出来画图”，而是建立一个稳定的明文 canonical dataset：读取层只处理 `subtask-c1-*` 和 `.SRFlux`，物理后处理层只处理频率积分和 source/transform 汇总，可视化层只处理 masked sampled field，从而避免目录结构耦合、频率拼合错误和插值造成的虚假 Poynting vector 分布。
