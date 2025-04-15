## To Run App Through Docker
0. Check that no virtual environment (venv) of any kind is running, you want to see the name of your computer as the first thing in the terminal 
1. docker-compose build --no-cache 
2. docker-compose up -d 
3. Paste "http://localhost:8000/" into browser and upload file to test 
4. If met with errors, paste "docker-compose logs frontend" into terminal and evaluate output 
5. To run code in Docker temporary environment, paste "docker exec -it scrna-frontend bash" into terminal and run "shiny run UI/anikaapp.py"


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

### Not needed anymore 
# for xg boost, libomp: first,
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" # takes about five mins
# add to path (for M chip mac)
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)"
# then
brew install libomp