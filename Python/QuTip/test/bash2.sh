#!/bin/bash
echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC10Nk6.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC20Nk4.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC20Nk6.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC20Nk8.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC20Nk10.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC10Nk15.py
echo "=== Analysis Ends ==="
