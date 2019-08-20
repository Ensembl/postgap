apt install -y git gcc make g++ curl autoconf pkg-config libz-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libgsl-dev libsqlite3-dev libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6 libssl-dev virtual-env

wget https://repo.anaconda.com/archive/Anaconda2-2019.07-Linux-x86_64.sh

sh Anaconda2-2019.07-Linux-x86_64.sh -b -p $PWD/anaconda
eval "$($PWD/anaconda/bin/conda shell.$0 hook)"
conda init
