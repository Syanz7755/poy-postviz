# SCUFF-EM Subtasks 通用拆解与目录生成规范

本文用于指导任意 geometry configuration 的 SCUFF-EM 任务拆解，不限定于
`c1`、`c4`、`d1` 或 refined mesh。后续接到新的 configuration 时，应先整理
配置清单，再按照本文生成可独立检查、独立提交和独立追踪结果的 subtask
directories。

当前 refined 案例的生成脚本位于：

```text
poy/scuffgeo/refined/generate_refined_tasks.py
```

该脚本是本文流程的一次具体实现；本文定义的流程和目录契约才是后续任务拆解的
通用参考。

## 1. 拆解模型

一个完整计算批次由两类变量组成：

1. Geometry configurations：网格、材料、相对位移及观测点设置的组合。
2. Frequency subsets：从完整材料频率数据中拆出的若干互不重叠子集。

每个 configuration 与每个 frequency subset 配对，形成一个独立 subtask：

```text
subtasks = configurations × frequency_subsets
```

若有 `C` 个 configurations、频率被拆为 `K` 个 subsets，则应生成：

```text
subtask_count = C × K
```

例如，3 个 configurations 与 4 个频率子集应生成 12 个 subtasks。

拆分后的每个 subtask 必须能够单独运行，不依赖其他 subtask 目录中的文件。

## 2. 开始前需要明确的信息

为每个新批次先填写以下批次参数：

| 参数 | 含义 | 示例 |
| --- | --- | --- |
| `batch_dir` | 本批次输出目录 | `poy/scuffgeo/refined` |
| `material_source` | 原始材料数据文件 | `poy/silica-27p.dat` |
| `selection_rule` | 是否对原始频率点降采样 | 从第 1 行开始隔行选取 |
| `subset_count` | 频率子集数量 | `4` |
| `configurations` | 本批次配置清单 | `c1`, `c4`, `d1` |
| `job_prefix` | Slurm job name 前缀 | `x` |
| `mesh_check_label` | 网格检查输出标签 | `meshcheck` |

然后为每个 configuration 填写一行配置清单：

| 字段 | 必填 | 说明 |
| --- | --- | --- |
| `config_id` | 是 | 简短且唯一的配置 ID，用于目录名和 job name |
| `lower_mesh` | 是 | Lower object 的源网格文件 |
| `upper_mesh` | 是 | Upper object 的源网格文件 |
| `displacement` | 是 | Upper object 的三维位移 `(dx, dy, dz)` |
| `evpoints_source` | 是 | 与该几何及位移匹配的 evpoints 文件 |
| `lower_material` | 是 | Lower 使用的完整材料插值表；当前 ordinary 任务统一为 `material.dat` |
| `upper_material` | 是 | Upper 使用的完整材料插值表；通常与 Lower 相同 |
| `notes` | 否 | 几何参数、网格精度、用途等说明 |

推荐先形成如下清单，再写或修改生成脚本：

| `config_id` | Lower mesh | Upper mesh | Displacement | evpoints |
| --- | --- | --- | --- | --- |
| `<id-a>` | `<path/to/lower>` | `<path/to/upper>` | `(0, 1.0, 0)` | `<path/to/evp>` |
| `<id-b>` | `<path/to/lower>` | `<path/to/upper>` | `(0, 4.0, 0)` | `<path/to/evp>` |

不要只根据文件名猜测配置关系。生成前应显式核对网格、位移和 evpoints
三者是否对应。

## 3. 命名规范

### 3.1 Configuration ID

`config_id` 应满足：

- 在同一批次中唯一。
- 只使用小写英文字母、数字、连字符或下划线。
- 足够简短，但能区分 geometry 和 displacement。
- 不包含 frequency subset 编号。

例如 `c1` 可以表示 C geometry、位移 1.0；若这种缩写可能产生歧义，可使用
更明确的 `c-dy1`。

### 3.2 Subtask directory

统一使用：

```text
subtask-<config_id>-<subset_index>
```

其中 `subset_index` 从 1 开始。示例：

```text
subtask-c1-1
subtask-c1-2
subtask-c4-1
subtask-d1-4
```

### 3.3 Slurm job name

job name 应短且唯一，推荐：

```text
<job_prefix><config_id><subset_index>
```

例如 `xc11` 表示 refined 批次中的 `c1` configuration、第 1 个频率子集。
如果多个批次会同时提交，应为不同批次使用不同 `job_prefix`，避免队列中难以区分。

## 4. 频率数据拆分

### 4.1 可选的降采样

若需要降低计算成本，先按任务要求从原始材料数据中选择频率点。选择时必须：

- 保留原始顺序。
- 记录明确且可复现的选择规则。
- 不修改每个数据点内部的列和值。
- 将完整的降采样结果另存为一个文件。

例如，从第 1 行开始每隔一行取一行：

```python
selected_lines = source_lines[0::2]
```

当前 refined 案例因此从 `silica-27p.dat` 得到 `silica-14p.dat`。

### 4.2 拆分为 subsets

设降采样后共有 `N` 个频率点，需要拆成 `K` 个 subsets。拆分结果必须满足：

```text
subset_i ∩ subset_j = ∅        对任意 i != j
subset_1 ∪ ... ∪ subset_K = selected_frequency_points
```

同时应满足：

- 每个频率点只出现一次。
- 所有子集合并后与降采样文件完全一致。
- 各 subset 尽量大小接近，以平衡计算时间。
- subset 内保持原始频率顺序。
- 不生成空 subset。

通用的均衡拆分方法如下：

```python
base, remainder = divmod(N, K)
sizes = [base + (i < remainder) for i in range(K)]
```

若任务另有“每个 subset 必须包含 2 至 3 个点”等限制，应先检查 `N` 和 `K`
是否能同时满足限制；无法满足时不得静默丢点，应调整 subset 数量或向任务定义者确认。

当前 refined 案例使用 14 个点和 4 个 subsets，实际大小为
`[3, 4, 3, 4]`。这是当前脚本的具体选择，不是其他批次必须照搬的固定比例。

### 4.3 生成 omegas

每个 frequency subset 对应一个 `omegas` 文件。必须区分以下两种频率单位：

- `material.dat` 第一列始终为物理角频率，单位是 `rad/s`。
- `scuff-neq --omegafile` 使用 SCUFF application omega 单位。输入物理角频率数值较大时，
  除以 `3e14`；已经是 SCUFF 单位的小数值不重复缩放。

当前生成器采用以下规则：

```python
SCUFF_OMEGA_SCALE = 3.0e14
SCUFF_OMEGA_LIMIT = 20.0

def format_scuff_omega(omega):
    value = omega / SCUFF_OMEGA_SCALE if abs(omega) >= SCUFF_OMEGA_LIMIT else omega
    if abs(value) >= SCUFF_OMEGA_LIMIT:
        raise ValueError("SCUFF omega must satisfy abs(omega) < 20")
    return f"{value:.12g}"

omega_lines = [
    format_scuff_omega(float(line.split()[0]))
    for line in frequency_subset
]
```

必须检查：

- `omegas` 行数等于当前 frequency subset 的行数。
- `omegas` 第 `i` 行与 frequency subset 第 `i` 行是同一物理频率，只是单位可能已缩放。
- 每个写入值均满足 `abs(omega) < 20`。
- 没有表头、空行或其他非频率内容被误写入。
- 不要缩放 `material.dat` 第一列；SCUFF 材料插值表仍使用 `rad/s`。

## 5. 生成每个 subtask

对每个 `(configuration, frequency_subset)` 组合执行以下步骤：

1. 创建 `subtask-<config_id>-<subset_index>`。
2. 创建其 `mshs` 目录。
3. 创建其 `logs` 目录。
4. 复制 Lower 和 Upper 网格到 `mshs`。
5. 准备 SCUFF-EM 实际引用的 `.msh` 文件。
6. 将完整 SCUFF-EM 材料插值表复制为 `material.dat`。
7. 从当前 frequency subset 生成并缩放 `omegas`。
8. 将 configuration 对应的 evpoints 复制为 `evpoints`。
9. 用 configuration 参数生成 `task.scuffgeo`。
10. 用 configuration ID 和 subset index 生成 `run.sbatch`。

推荐的目录契约为：

```text
subtask-<config_id>-<subset_index>/
|-- evpoints
|-- omegas
|-- run.sbatch
|-- material.dat
|-- task.scuffgeo
|-- logs/
`-- mshs/
    |-- lower.msh
    `-- upper.msh
```

源网格若为 `.nas`，但稍后会转为同名 `.msh`，可以同时保留：

```text
mshs/
|-- <lower-name>.nas
|-- <lower-name>.msh
|-- <upper-name>.nas
`-- <upper-name>.msh
```

注意：仅修改扩展名并不等同于格式转换。若 SCUFF-EM 需要真实的 Gmsh 格式，
必须在运行计算前完成实际转换，并用 `scuff-analyze` 验证。

## 6. `task.scuffgeo` 通用模板

```text
OBJECT Lower
  MESHFILE mshs/<lower-mesh>.msh
  MATERIAL FILE_material.dat
ENDOBJECT

OBJECT Upper
  MESHFILE mshs/<upper-mesh>.msh
  MATERIAL FILE_material.dat
  DISPLACED <dx> <dy> <dz>
ENDOBJECT
```

生成时必须替换以下占位符：

| 占位符 | 来源 |
| --- | --- |
| `<lower-mesh>` | configuration 的 `lower_mesh` |
| `<upper-mesh>` | configuration 的 `upper_mesh` |
| `<dx> <dy> <dz>` | configuration 的 `displacement` |

若 Lower 和 Upper 使用不同材料，应分别复制完整材料表并使用不同文件名，不要继续共用
`material.dat`。

## 7. `run.sbatch` 通用模板

```bash
#!/bin/bash
#SBATCH --job-name=<job-name>
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=logs/<subtask-name>-%j.out
#SBATCH --error=logs/<subtask-name>-%j.err

mkdir -p logs

module load scuff-em

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false

scuff-neq \
  --geometry task.scuffgeo \
  --epfile evpoints \
  --omegafile omegas \
  --emtpft \
  --sourceobject Lower \
  --destobject Upper
```

模板中的文件路径必须与 subtask 的实际目录结构一致。若选择把 evpoints 放进
子目录，例如 `evpoints/<name>.dat`，则 `--epfile` 必须同步修改，不能混用两种约定。
`<subtask-name>` 应替换为完整目录名，例如 `subtask-cd1-1`；`%j` 由 Slurm 替换为
job ID，因此日志名形如 `logs/subtask-cd1-1-59322525.out`。

partition、CPU 数、module 名称和 `scuff-neq` 参数属于计算环境设置。迁移到其他
集群时应单独检查，不应由 configuration ID 推导。

## 8. 批量检查与提交

在批次根目录放置统一的网格检查脚本：

```bash
#!/bin/bash
set -euo pipefail

for dir in subtask-*; do
  [ -d "$dir" ] || continue
  (
    cd "$dir"
    scuff-analyze \
      --geometry task.scuffgeo \
      --writegmshfiles meshcheck
  )
done
```

运行：

```bash
bash analyze_all_scuffgeo.sh
```

检查通过后，再使用批量提交脚本：

```bash
#!/bin/bash
set -euo pipefail

for dir in subtask-*; do
  [ -d "$dir" ] || continue
  (
    cd "$dir"
    sbatch run.sbatch
  )
done
```

运行：

```bash
bash submit_all_slurm.sh
```

不要在尚未完成输入验证和网格检查时直接批量提交。

## 9. 自动验证要求

生成脚本完成后，应至少自动检查以下项目。

### 9.1 批次级检查

- 实际 subtask 数等于 `configuration_count × subset_count`。
- 所有 `config_id` 唯一。
- 所有源网格、材料文件和 evpoints 文件存在。
- 所有 subsets 无重叠、无遗漏，合并后等于降采样结果。
- 批量检查和批量提交脚本存在。

### 9.2 Subtask 级检查

- `task.scuffgeo`、`run.sbatch`、`evpoints`、`omegas` 和 `material.dat` 均存在。
- `logs` 和 `mshs` 目录均存在。
- `mshs` 中 Lower 和 Upper 网格均存在。
- `task.scuffgeo` 引用的文件都位于当前 subtask 内。
- `DISPLACED` 与 configuration 清单一致。
- `--epfile` 与实际 evpoints 路径一致。
- `--omegafile` 与实际 omegas 路径一致。
- `omegas` 与 frequency subset 行数相同、顺序一致，且所有值的绝对值小于 20。
- `material.dat` 第一列保留物理角频率 `rad/s`，没有执行 `3e14` 缩放。
- stdout/stderr 路径均为 `logs/<subtask-name>-%j.<suffix>`。
- Slurm job name 在本批次内唯一。

### 9.3 计算环境检查

- `scuff-analyze` 能读取每个 `task.scuffgeo`。
- 网格转换真实有效，不只是改扩展名。
- Slurm partition 和资源参数在目标集群上可用。
- `module load scuff-em` 能成功加载所需版本。

## 10. 推荐的生成器结构

为了让一个脚本适用于后续不同 configurations，建议将“变化的数据”和“固定的
生成逻辑”分离。配置可集中写成：

```python
CONFIGS = {
    "<config-id>": {
        "lower_mesh": ROOT / "in" / "<lower-file>",
        "upper_mesh": ROOT / "in" / "<upper-file>",
        "displacement": ("<dx>", "<dy>", "<dz>"),
        "evpoints": ROOT / "evpoints" / "<evpoints-file>",
    },
}
```

固定生成逻辑只遍历配置和频率子集：

```python
for config_id, config in CONFIGS.items():
    for subset_index, subset in enumerate(subsets, start=1):
        create_subtask(config_id, config, subset_index, subset)
```

新增 configuration 时，原则上只增加一条配置记录，不复制整段目录生成代码。
当 configurations 较多时，可进一步将清单移到 JSON、YAML 或 CSV 文件，使生成器
不需要因新增配置而修改核心逻辑。

## 11. 新 configuration 执行清单

后续处理任意新 configuration 时，按以下顺序执行：

- [ ] 为 configuration 选定唯一 `config_id`。
- [ ] 确认 Lower 和 Upper 网格路径。
- [ ] 确认三维 displacement。
- [ ] 确认对应 evpoints，不依赖文件名猜测。
- [ ] 确认材料文件和频率选择规则。
- [ ] 确认 subset 数量及每个 subset 的大小约束。
- [ ] 生成 configuration 与 frequency subsets 的笛卡尔积。
- [ ] 验证所有 subtask 的目录契约和文件引用。
- [ ] 运行全部 `scuff-analyze` 检查。
- [ ] 确认检查通过后再批量 `sbatch`。
- [ ] 保存本批次配置清单和生成脚本，以便复现。

## 12. 当前 refined 案例

当前 `poy/scuffgeo/refined` 只是上述通用流程的一个实例：

| `config_id` | 网格 | Displacement | evpoints |
| --- | --- | --- | --- |
| `c1` | C refined mesh | `(0, 1.0, 0)` | `c_d1.0.txt` |
| `c4` | C refined mesh | `(0, 4.0, 0)` | `c_d4.0.txt` |
| `d1` | D refined mesh | `(0, 1.0, 0)` | `d_d1.0.txt` |

材料数据由 `silica-27p.dat` 隔行选取得到 14 点，再拆成 4 个 subsets，因此生成
12 个 subtask directories。后续其他 configuration 可以使用不同网格、位移、
evpoints、材料数据、降采样规则和 subset 数量，但应保持本文定义的拆解关系、
目录自包含原则和验证流程。

## 13. 当前 ordinary `cd1` / `cd4` 生成方式

生成器位于：

```text
poy/scuffgeo/ordinary/generate_ordinary_tasks.py
```

只重新生成 `cd1` 和 `cd4`：

```powershell
python scuffgeo\ordinary\generate_ordinary_tasks.py --configs cd1 cd4
```

当前材料源为：

```text
D:\Researches\wqs_comps\aligned_SiO2-Franta-300C.txt
```

源文件三列依次为 `omega(rad/s)`、`Re(epsilon)`、`Im(epsilon)`。生成器将其转换为
SCUFF-EM 的两字段格式 `omega complex-epsilon`，并把完整 2502 点表写入每个
subtask 的 `material.dat`。计算频率仍从 `silica-27p.dat` 选择并拆成 9 个 subsets；
写入 `omegas` 和批次级 `omegas-sub*.dat` 时按第 4.3 节转换为 SCUFF application
omega 单位。

生成器会报告 calculation omega 是否超出 `material.dat` 的插值频率范围。当前
27 点频率表有 8 点低于新材料表的最小频率；这不会被静默删除或替换，正式计算前
应根据物理需求决定扩展材料表或调整 calculation omega。
