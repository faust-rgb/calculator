import subprocess

tests = [
    ("高次有理分式 1", "integral(1 / (x^4 + 1))"),
    ("高次有理分式 2", "integral(x / (x^4 + 2 * x^2 + 2))"),
    ("高次有理分式 3", "integral((x^2 + 1) / (x^4 + 1))"),
    ("双二次有理分式", "integral(1 / ((x^2 + 1) * (x^2 + 4)))"),
    ("嵌套指数微分塔", "integral(exp(exp(x)) * exp(x))"),
    ("嵌套对数微分塔", "integral(ln(ln(x)) / x)"),
    ("复合指数代换", "integral(x * exp(x^2))"),
    ("多项式-对数乘积", "integral(x^2 * ln(x))"),
    ("指数-正弦复合积分", "integral(sin(x) * exp(x))"),
    ("多项式-余弦复合积分", "integral(x * cos(x))"),
    ("正切积分", "integral(tan(x))"),
    ("正弦奇次幂积分", "integral(sin(x)^3)"),
    ("根式代数扩展 1", "integral(x * sqrt(x + 1))"),
    ("根式代数扩展 2", "integral(x / sqrt(1 - x^2))"),
    ("欧拉换元代数扩展", "integral(1 / (x * sqrt(x^2 + 1)))"),
    ("严格非初等/指数积分", "integral(exp(x) / x)"),
    ("严格非初等/高斯误差函数", "integral(exp(-x^2))"),
    ("严格非初等/正弦积分", "integral(sin(x) / x)"),
    ("严格非初等/余弦积分", "integral(cos(x) / x)"),
    ("严格非初等/对数积分", "integral(1 / ln(x))"),
    ("代数曲线非初等/椭圆积分", "integral(1 / sqrt(x^3 + 1))"),
    ("代数曲线非初等/超椭圆积分", "integral(1 / sqrt(x^5 - x + 1))")
]

print("-" * 125)
print(f"| {'类别':<22} | {'输入表达式':<35} | {'CLI 计算结果':<60} |")
print("-" * 125)

for category, expr in tests:
    input_str = f":symbolic on\n{expr}\nexit\n"
    proc = subprocess.run(["./bin/calculator"], input=input_str, capture_output=True, text=True)
    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip() and not line.startswith("Command Line") and not line.startswith("Enter") and not line.startswith(">")]
    res = lines[-1] if lines else "Error"
    if res.startswith(">"):
        res = res[1:].strip()
    print(f"| {category:<20} | {expr:<35} | {res:<60} |")

print("-" * 125)
