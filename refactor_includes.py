import os
import re

src_root = "/home/lenovo/code/calculator/src"
test_root = "/home/lenovo/code/calculator/test"

def get_header_map(root, rel_base):
    header_map = {} # filename -> list of relative paths
    full_path_set = set() # set of all relative paths
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if f.endswith('.h'):
                rel_path = os.path.relpath(os.path.join(dirpath, f), rel_base)
                if f not in header_map:
                    header_map[f] = []
                header_map[f].append(rel_path)
                full_path_set.add(rel_path)
    return header_map, full_path_set

src_headers, src_full_paths = get_header_map(src_root, src_root)
test_headers, test_full_paths = get_header_map(test_root, test_root)

# Layer definitions (lower to higher)
layer_hierarchy = [
    {"math"},
    {"types", "matrix"},
    {"io"},
    {"parser"},
    {"execution"},
    {"core"},
    {"plot", "statistics", "dsp", "symbolic", "polynomial", "time"},
    {"analysis"},
    {"module"},
    {"app"}
]

layer_to_rank = {}
for i, group in enumerate(layer_hierarchy):
    for layer in group:
        layer_to_rank[layer] = i

def get_layer(path):
    parts = path.split(os.sep)
    if not parts:
        return None
    if parts[0] in layer_to_rank:
        return parts[0]
    return None

violations = []

def process_file(filepath, is_src=True):
    with open(filepath, 'r') as f:
        content = f.read()

    rel_filepath = os.path.relpath(filepath, src_root if is_src else test_root)
    file_layer = get_layer(rel_filepath) if is_src else "test"

    def replace_include(match):
        inc_path = match.group(1)
        
        # Check if it's already a full path in src
        if inc_path in src_full_paths:
            actual_inc = inc_path
        elif os.path.basename(inc_path) in src_headers:
            # It's a short path, or a path that needs fixing
            filename = os.path.basename(inc_path)
            candidates = src_headers[filename]
            if len(candidates) == 1:
                actual_inc = candidates[0]
            else:
                file_dir = os.path.dirname(rel_filepath)
                same_dir_candidates = [c for c in candidates if os.path.dirname(c) == file_dir]
                if same_dir_candidates:
                    actual_inc = same_dir_candidates[0]
                else:
                    actual_inc = candidates[0]
        elif inc_path in test_full_paths:
            actual_inc = inc_path
            return f'#include "{actual_inc}"' # No layering check for test files in this way
        elif os.path.basename(inc_path) in test_headers:
            actual_inc = test_headers[os.path.basename(inc_path)][0]
            return f'#include "{actual_inc}"'
        else:
            return match.group(0)

        # Check for layering violation
        if is_src and file_layer:
            inc_layer = get_layer(actual_inc)
            if inc_layer and file_layer in layer_to_rank and inc_layer in layer_to_rank:
                if layer_to_rank[inc_layer] > layer_to_rank[file_layer]:
                    violations.append(f"Violation in {rel_filepath} (layer {file_layer}): includes {actual_inc} (layer {inc_layer})")
        
        return f'#include "{actual_inc}"'

    new_content = re.sub(r'#include\s+"([^"]+)"', replace_include, content)
    
    if new_content != content:
        with open(filepath, 'w') as f:
            f.write(new_content)
        return True
    return False

files_to_process = []
for root, _, filenames in os.walk(src_root):
    for f in filenames:
        if f.endswith(('.cpp', '.h')):
            files_to_process.append((os.path.join(root, f), True))

for root, _, filenames in os.walk(test_root):
    for f in filenames:
        if f.endswith(('.cpp', '.h')):
            files_to_process.append((os.path.join(root, f), False))

for filepath, is_src in files_to_process:
    process_file(filepath, is_src)

# Report violations uniquely
unique_violations = sorted(list(set(violations)))
print(f"Found {len(unique_violations)} unique layering violations:")
for v in unique_violations:
    print(v)
