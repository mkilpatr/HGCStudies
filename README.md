### HGCStudies

Hexacontroller Test Stand

First Test stand IP: 128.111.19.46

Second Test stand IP: 128.111.19.45

Load FW to FPGA:
```
ssh root@baryon.physics.ucsb.edu # login from hgcal PC
cd mylittledt
make distclean
make trophy
```
start processes for running
Locally run Fast Control:
```
source env.sh
cd /root/hexactrl-sw
./bin/zmq-server &
```

Locally run Slow Control:
```
cd /root/hexactrl-sw/zmq_i2c
python3 zmq_server.py &
```

Need two separate terminals logged into hgcal pc

In first terminal:
```
cd hexactrl-sw/
source env.sh
./bin/zmq-client
```

In second terminal:
```
cd hexactrl-sw/hexactrl-script/
source env.sh
python3 full_test.py -d <outputdir> -i <test stand IP> -f configs/init<HD or LD depending on module>.yaml
```
