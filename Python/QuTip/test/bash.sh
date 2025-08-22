#!/bin/bash
echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC10Nk2.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC10Nk4.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC10Nk8.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC10Nk10.py
echo "=== Analysis Ends ==="

echo "=== HEOM Performance Analysis ==="
/usr/bin/time -v python NC20Nk4.py
echo "=== Analysis Ends ==="