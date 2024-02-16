Bootstrap: docker
From: python:3.12

%files
./kraken2ref/ /home/kraken2ref/kraken2ref
./.git /home/kraken2ref/.git
./pyproject.toml /home/kraken2ref/
./setup.cfg  /home/kraken2ref/
./requirements.txt  /home/kraken2ref/

%post
cd /home/kraken2ref/ && pip install .
