# Install GCC7 For Ubuntu 16.04
#
# Run with: source installGCC7.sh

sudo apt-get install -y software-properties-common
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install g++-7 -y
