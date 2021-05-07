# Launch ray cluster
#
# Modify lines below according to your cluster
# 
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8


# Head
HEAD_IP=192.168.10.201

# Nodes
IP_ARRAY=(192.168.10.200 192.168.10.202 192.168.10.203)

# Conda virtual environment on each node
CONDA_ENV=graniitti

# Paths
LOCAL_DIR=${HOME}/cernbox/graniitti/python
REMOTE_DIR=${HOME}/cernbox/graniitti
CONDA_PATH=${HOME}/anaconda3

# ------------------------------------------------------------------------

port=6379
address="$HEAD_IP:$port"
echo $address


conda activate $CONDA_ENV
ray stop --force
sleep 2
ray stop --force
ray start --head --port=$port

# Launch remote nodes

# pip install -r $REMOTE_DIR/requirements.txt; \

cmd=". ~/.bashrc; source $CONDA_PATH/etc/profile.d/conda.sh; conda activate $CONDA_ENV; \
chmod +x $REMOTE_DIR/bin/*; \
ray stop --force; sleep 2; ray stop --force; ray start --address=$address"


for IP in "${IP_ARRAY[@]}"
do
	echo "connecting remote node: ${IP}"
	
	# Sync the files (rsync will copy only files which changed)
	rsync -avz -e 'ssh' ${LOCAL_DIR} user@${IP}:${REMOTE_DIR}

	# Start ray
	ssh ${IP} ${cmd}
done

ray status
