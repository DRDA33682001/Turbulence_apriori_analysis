import re

with open('REPORT.md', 'r') as f:
    text = f.read()

matches = re.findall(r'\*\*(?:Figure|Table):\*\*.*?(?:\n\n\$|\n\$).*?\}', text, re.DOTALL)
for i, m in enumerate(matches):
    print(f"--- Match {i} ---")
    print(m)

