### HGCStudies ###

## Hexacontroller Test Stand ##

First Test stand IP: 128.111.19.46

Second Test stand IP: 128.111.19.45

# Load FW to FPGA #
```
ssh root@baryon.physics.ucsb.edu # login from hgcal PC
newData # needed to set the date of the PC to something that's not 1970
cd mylittledt
make distclean
make trophy
```
# Start processes for running #
```
cd /root/hexactrl-sw
python3 address_table/buildXML.py
python3 webserver.py & # not necessary if you plan to run within the terminal
```

# Locally run Fast Control #
```
source env.sh
cd /root/hexactrl-sw
./bin/zmq-server &
```

# Locally run Slow Control #
```
cd /root/hexactrl-sw/zmq_i2c
python3 zmq_server.py &
```

## Need two separate terminals logged into hgcal pc ##

# In first terminal #
```
cd hexactrl-sw/
source env.sh
./bin/zmq-client
```

# In second terminal #
```
cd hexactrl-sw/hexactrl-script/
source env.sh
python3 full_test.py -d <outputdir> -o /data/HexaController/data/ -i <test stand IP> -f configs/<initHD.yaml or initLD.yaml depending on module>
```
