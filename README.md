# scRNA-seq-Automation

pip install -r requirements.txt

python src/main.py

shiny run UI/anikaapp.py

docker-compose build --no-cache
after this, you need to run "docker-compose up -d" to actually start the docker
and when leaving VSCode, run "docker-compose down"  # To run in background

git add <file>        # After resolving
git commit -m "<message>" # Finish merge
git push origin main   # Push if needed

git status # to see what files you have made changes to