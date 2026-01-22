#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi

#!/usr/bin/env bash
# ------------------------------------------------------------
# Checks on a directory containing *.fastq.gz files
#   • Directory must exist
#   • Directory must not be empty
#   • No more than ONE empty fastq.gz file is allowed
#   • Every fastq.gz file must be at least 80 bytes
# ------------------------------------------------------------

# Assign the first argument to the DIRECTORY variable
DIRECTORY="$1"

# ---------- 1. Directory existence ----------
if [[ -z "$DIRECTORY" ]]; then
    echo "Error: No directory argument supplied."
    exit 1
fi

if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# ---------- 2. Directory non‑empty ----------
if [[ -z "$(ls -A "$DIRECTORY")" ]]; then
    echo "Error: Directory '$DIRECTORY' is empty."
    exit 1
fi

# ---------- 3. Find empty *.fastq.gz files ----------
empty_files=$(find "$DIRECTORY" -type f -name "*.fastq.gz" -empty)

# Count how many empty fastq.gz files were found
empty_count=$(echo "$empty_files" | grep -c '^' || true)   # 0 if none

if (( empty_count > 1 )); then
    echo "Error: More than one empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
    exit 1
fi

# ---------- 4. Check file size (>= 80 bytes) ----------
# Loop over all fastq.gz files and verify size
while IFS= read -r -d '' file; do
    size=$(stat -c%s "$file")   # size in bytes
    if (( size < 80 )); then
        echo "Error: File '$file' is smaller than 80 bytes (size=$size)."
        exit 1
    fi
done < <(find "$DIRECTORY" -type f -name "*.fastq.gz" -print0)

# ---------- 5. Success message ----------
if [[ -n "$empty_files" ]]; then
    echo "Warning: One empty *.fastq.gz file found in '$DIRECTORY':"
    echo "$empty_files"
else
    echo "All *.fastq.gz files in '$DIRECTORY' are present and ≥ 80 bytes."
fi
