import os
import shutil
import time

src_base = r"G:\내 드라이브\ITLAS\results"
dest_dir = r"G:\내 드라이브\ITLAS\g_paper-IT-v16_26f23-\g_paper-IT-v16_image-26f25"

suffix_map = {
    "nb01_qc": "_qd",
    "nb02_composition": "_com",
    "nb04_global": "_glo",
    "nb05_lineage": "_lin",
    "nb06_subcluster": "_sub",
    "nb08_correlation": "_cor",
    "nb09_figures": "_fig",
    "nb10_supplementary": "_sup",
    "task3_donor_validation": "_don",
}

os.makedirs(dest_dir, exist_ok=True)
image_exts = {".png", ".jpg", ".jpeg", ".svg", ".pdf", ".tiff", ".tif"}
copied_count = 0
skipped_count = 0

for root, dirs, files in os.walk(src_base):
    for file in files:
        ext = os.path.splitext(file)[1].lower()
        if ext in image_exts:
            rel_path = os.path.relpath(root, src_base)
            top_level_dir = rel_path.split(os.sep)[0]

            suffix = suffix_map.get(top_level_dir, "")

            name, ex = os.path.splitext(file)
            new_name = f"{name}{suffix}{ex}"
            src_path = os.path.join(root, file)
            dest_path = os.path.join(dest_dir, new_name)

            if os.path.exists(dest_path):
                try:
                    if os.path.getsize(src_path) == os.path.getsize(dest_path):
                        skipped_count += 1
                        continue
                except OSError:
                    pass

            max_retries = 3
            for attempt in range(max_retries):
                try:
                    shutil.copy2(src_path, dest_path)
                    copied_count += 1
                    print(f"Copied {file} -> {new_name}")
                    break
                except Exception as e:
                    print(
                        f"Error copying {file}: {e}. Retrying {attempt + 1}/{max_retries}"
                    )
                    time.sleep(2)
            else:
                print(f"Failed to copy {file} after {max_retries} retries.")

print(f"\nTotal copied: {copied_count} images. Skipped: {skipped_count}")
