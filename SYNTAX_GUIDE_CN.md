# 计算器脚本语法指南（中文版）

本文档是 `.calc` 脚本语法的完整中文参考。脚本采用 Python 风格的语法：`#` 注释、冒号起始的代码块、缩进定义作用域、语句末尾无需分号。

## 运行脚本

从项目根目录运行：

```bash
./bin/calculator script.calc
```

在交互式提示符中运行：

```text
> :run script.calc
```

脚本文件应使用 `.calc` 扩展名。

## 文件结构

- 使用 `#` 进行行注释
- 每行一条语句
- 语句末尾不需要 `;`
- `if`、`elif`、`else`、`while`、`for`、`def` 代码块以 `:` 开始
- 使用一致的缩进定义代码块体，推荐使用空格
- 矩阵字面量内部仍使用 `;` 分隔行，例如 `[1, 2; 3, 4]`

```calc
# 这是一个注释
x = 10
print("x = ", x)
```

## 值与变量

脚本支持数字、字符串、标量变量、表达式风格的自定义函数和矩阵。

```calc
count = 3
label = "ready"
m = [1, 2; 3, 4]
n = mat(2, 2, 5, 6, 7, 8)
```

变量名区分大小写。`pi`、`e` 等常数、内置函数名和活动函数名不应被重新定义为普通变量。

当变量需要在代码块之后被读取时，应在代码块之前初始化，并在代码块内部赋值。

### 关键字和保留名

脚本关键字不能用作变量名：

```
def, fn, if, elif, else, while, for, in, return, break, continue, pass
```

保留函数名不能用于用户定义的函数。完整列表见 `KEYWORDS_REFERENCE.md`。常见示例：

```
sin, cos, tan, exp, ln, sqrt, abs, det, trace, inv, transpose, diff, integral
```

变量命名规则：

- 首字符必须是字母或下划线
- 后续字符可以是字母、数字或下划线
- 不能是关键字或保留函数名

## 表达式

计算器表达式使用与 REPL 相同的运算符和内置函数。

```calc
total = 1 + 2 * 3
power = 2 ^ 10
distance = sqrt(3 ^ 2 + 4 ^ 2)
ok = abs(sin(pi / 2) - 1) < 0.000001
```

在脚本条件中，使用嵌套的 `if` 语句代替 `&&` 或 `||`。

## 控制流

```calc
if score >= 90:
    grade = "A"
elif score >= 80:
    grade = "B"
else:
    grade = "other"

i = 1
sum = 0
while i <= 10:
    sum = sum + i
    i = i + 1

odd_sum = 0
for j in range(0, 10):
    if j == 7:
        break
    if mod(j, 2) == 0:
        continue
    odd_sum = odd_sum + j
```

支持 `range(stop)`、`range(start, stop)` 和 `range(start, stop, step)`。循环在 `stop` 之前停止，与 Python 约定一致。

## 脚本函数

使用 `def` 定义脚本函数。当函数需要返回值时，必须使用 `return`。

```calc
def fact_script(n):
    if n <= 1:
        return 1
    return n * fact_script(n - 1)

fact_script(6)
```

避免使用已作为内置函数名的名称定义脚本函数。

## 表达式风格函数

表达式风格函数可在脚本中定义，然后用于分析、多项式和微积分命令。

```calc
quad_fn(x) = x ^ 2 + 2 * x + 1
diff(quad_fn)
integral(quad_fn, 0, 1)
taylor(quad_fn, 0, 3)
```

## 输出

`print(...)` 接受脚本值、字符串和格式化矩阵。符号命令如 `diff(x ^ 2)` 或 `laplace(step(t), t, s)` 通常应作为独立语句编写，而不是嵌套在 `print(...)` 内部。

```calc
print("det = ", det([1, 2; 3, 4]))
laplace(step(t), t, s)
```

## 性能优化

脚本引擎包含多项性能优化：

### 表达式缓存

循环体内的表达式会被自动缓存。首次执行时编译为 AST，后续迭代直接复用：

```calc
# 表达式 (a * b + c) * (i + 1) - i 只编译一次
result = 0
a = 2
b = 3
c = 5
for i in range(0, 10000):
    result = result + (a * b + c) * (i + 1) - i
```

### 变量查找优化

局部变量使用平坦数组存储，查找效率高于传统的红黑树实现：

- 函数参数和局部变量存储在连续内存中
- 变量查找从内层作用域向外层搜索
- 编译期可将变量绑定到固定槽位索引

### 栈帧优化

函数调用使用预分配的栈帧机制，避免每次调用都动态分配内存：

```calc
# 递归函数调用效率更高
def my_fib(n):
    if n <= 1:
        return n
    return my_fib(n - 1) + my_fib(n - 2)

my_fib(20)  # 快速执行
```

## 兼容性说明

解析器内部仍接受旧的花括号加分号形式，以便加载包含脚本函数的已保存状态文件。新脚本和示例应使用上述 Python 风格形式。

## 功能覆盖示例

使用 `test/script/comprehensive_validation.calc` 作为广泛的回归测试脚本。它覆盖：

- 注释、赋值、字符串、算术和比较
- `if/elif/else`、嵌套条件、`while`、`for range`、`break` 和 `continue`
- 递归 `def` 函数和表达式风格函数
- 数值数学辅助函数和常数
- 矩阵、行列式、迹、秩、逆、求解、转置和 rref
- 多项式向量辅助函数
- 导数、积分、极限、泰勒、求根、简化、临界点、Pade、拉普拉斯、傅里叶和 Z 变换
- 程序员风格的位辅助函数

预期最终输出：

```text
11
```

## 示例脚本

```calc
# 综合示例脚本

print("== 计算器脚本示例 ==")

# 基本运算
a = 10
b = 20
c = a + b
print("a + b = ", c)

# 条件语句
if c == 30:
    print("条件成立")
else:
    print("条件不成立")

# while 循环
sum = 0
i = 1
while i <= 10:
    sum = sum + i
    i = i + 1
print("1+2+...+10 = ", sum)

# for 循环
for_sum = 0
for j in range(1, 11):
    for_sum = for_sum + j
print("for 循环求和 = ", for_sum)

# 函数定义
def square(x):
    return x * x

print("square(5) = ", square(5))

# 递归函数
def factorial(n):
    if n <= 1:
        return 1
    return n * factorial(n - 1)

print("5! = ", factorial(5))

# 矩阵操作
m = [1, 2; 3, 4]
print("矩阵 m = ", m)
print("行列式 = ", det(m))
print("逆矩阵 = ", inverse(m))

# 微积分
f(x) = x ^ 2 + 2 * x + 1
print("f(x) = ", f(x))
print("导数 = ", diff(f))
print("定积分 [0,1] = ", integral(f, 0, 1))

print("== 示例结束 ==")
"完成"
```
