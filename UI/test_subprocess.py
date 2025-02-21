import subprocess

# Test invoking Python via subprocess
result = subprocess.run(
    ['python', '-c', 'print("Hello from subprocess!")'],
    capture_output=True, text=True
)

print(f"stdout: {result.stdout}")
print(f"stderr: {result.stderr}")
