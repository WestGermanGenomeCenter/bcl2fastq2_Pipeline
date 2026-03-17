import csv
import argparse
import sys
import re
from collections import Counter

class SampleSheetValidator:
    def __init__(self, file_path):
        self.file_path = file_path
        self.sections = {}
        self.errors = []
        self.warnings = []

    def parse(self):
        """Parses the SampleSheet into sections."""
        current_section = None
        try:
            with open(self.file_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                for line in reader:
                    if not line or all(x == '' for x in line):
                        continue
                    
                    header_match = re.match(r'^\[(.*)\]', line[0])
                    if header_match:
                        current_section = header_match.group(1)
                        self.sections[current_section] = []
                    elif current_section:
                        self.sections[current_section].append(line)
        except Exception as e:
            self.errors.append(f"CRITICAL FILE ERROR: Could not read file. {str(e)}")

    def validate(self):
        """Performs BCL Convert specific validations."""
        data_key = 'BCLConvert_Data' if 'BCLConvert_Data' in self.sections else 'Data'
        
        if data_key not in self.sections:
            self.errors.append("MISSING REQUIRED SECTION: [BCLConvert_Data]")
            return

        data_rows = self.sections[data_key]
        if not data_rows:
            self.errors.append(f"SECTION [{data_key}] IS EMPTY.")
            return

        header = [h.strip() for h in data_rows[0]]
        rows = data_rows[1:]

        # Required columns for BCL Convert
        required_cols = ['Sample_ID', 'Index']
        for col in required_cols:
            if col not in header:
                self.errors.append(f"MISSING REQUIRED COLUMN in [{data_key}]: {col}")
        
        if self.errors: return 

        col_map = {col: i for i, col in enumerate(header)}
        
        # Trackers
        sample_ids = []
        lane_index_combos = {} # {lane: [index1+index2, ...]}

        for row_idx, row in enumerate(rows, start=2):
            if len(row) < len(header):
                self.errors.append(f"LINE {row_idx}: Row has fewer columns than header.")
                continue

            s_id = row[col_map['Sample_ID']].strip()
            idx1 = row[col_map['Index']].strip().upper()
            idx2 = row[col_map.get('Index2', 0)].strip().upper() if 'Index2' in col_map else ""
            lane = row[col_map.get('Lane', 0)].strip() if 'Lane' in col_map else "all"

            # 1. Validate Sample_ID
            if not re.match(r'^[a-zA-Z0-9_-]+$', s_id):
                self.errors.append(f"LINE {row_idx}: Invalid Sample_ID '{s_id}'. Use only alphanumeric, dashes, or underscores.")

            # 2. Validate Index Sequences
            if not re.match(r'^[ATCGN]*$', idx1):
                self.errors.append(f"LINE {row_idx}: Invalid characters in Index '{idx1}'.")
            if idx2 and not re.match(r'^[ATCGN]*$', idx2):
                self.errors.append(f"LINE {row_idx}: Invalid characters in Index2 '{idx2}'.")

            # 3. Validate Adapter Sequences (if present)
            for adapter_col in ['Adapter', 'AdapterRead2']:
                if adapter_col in col_map:
                    seq = row[col_map[adapter_col]].strip().upper()
                    if seq and not re.match(r'^[ATCG]*$', seq):
                        self.errors.append(f"LINE {row_idx}: Invalid DNA bases in {adapter_col} '{seq}'.")

            sample_ids.append(s_id)
            
            # Store index combos per lane
            if lane not in lane_index_combos:
                lane_index_combos[lane] = []
            lane_index_combos[lane].append(f"{idx1}+{idx2}")

        # 4. Check for Duplicate Sample IDs (Global)
        id_counts = Counter(sample_ids)
        for sid, count in id_counts.items():
            if count > 1:
                self.errors.append(f"DUPLICATE SAMPLE_ID: '{sid}' appears {count} times.")

        # 5. Check for Index Clashes (Per Lane)
        for lane, combos in lane_index_combos.items():
            idx_counts = Counter(combos)
            for combo, count in idx_counts.items():
                if count > 1:
                    lane_info = f" in Lane {lane}" if lane != "all" else ""
                    self.errors.append(f"INDEX COLLISION{lane_info}: Combination '{combo}' appears {count} times.")

    def report(self):
        """Prints results to console."""
        if self.errors:
            print(f"[FAIL] Validation Failed for {self.file_path}:")
            for err in self.errors:
                print(f"  - {err}")
            return 0 
        else:
            print(f"[SUCCESS] SampleSheet {self.file_path} is valid for BCL Convert.")
            if self.warnings:
                for warn in self.warnings:
                    print(f"  [WARN] {warn}")
            return 0

def main():
    parser = argparse.ArgumentParser(description="Illumina BCL Convert SampleSheet Validator")
    parser.add_argument("samplesheet", help="Path to the SampleSheet.csv")
    args = parser.parse_args()

    validator = SampleSheetValidator(args.samplesheet)
    validator.parse()
    validator.validate()
    sys.exit(validator.report())

if __name__ == "__main__":
    main()