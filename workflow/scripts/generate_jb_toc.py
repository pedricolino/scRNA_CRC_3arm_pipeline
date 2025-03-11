import os

notebooks_dir = "jupyter_book/chosen_branch"  # Change this to your notebooks folder
toc_file = "jupyter_book/_toc.yml"

# Header for _toc.yml
toc_content = """format: jb-book
root: intro  # Change this to your actual root file
chapters:
"""

# List all `.ipynb` files in the folder (sorted alphabetically)
notebooks = sorted([f for f in os.listdir(notebooks_dir) if f.endswith(".ipynb")])

# Add each notebook to the `_toc.yml` structure
for nb in notebooks:
    toc_content += f"  - file: {notebooks_dir.replace("jupyter_book/", "")}/{nb.replace('.ipynb', '')}\n"

# Write to _toc.yml
with open(toc_file, "w") as f:
    f.write(toc_content)

print(f"Generated {toc_file} with {len(notebooks)} notebooks.")