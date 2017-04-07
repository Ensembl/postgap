## BASIC INSTRUCTIONS
#####################

Start server on http://0.0.0.0:8080/ :

python server.py

Edit server.py to change host/port.



## ENDPOINTS
############

Currently there is only one endpoint /query that can be accessed in a few ways:

# By rsID
http://0.0.0.0:8080/query?rsID=rs11712165

# By rsID, chr and pos
http://0.0.0.0:8080/query?rsID=rs11712165&chr=1&pos=230845794

# By EFO (allows multiple)
http://0.0.0.0:8080/query?efos=EFO_0000408
http://0.0.0.0:8080/query?efos=EFO_0000408&efos=EFO_0000409

# By disease (allows multiple)
http://0.0.0.0:8080/query?diseases="breast%20cancer"



## DEBUGGING
############

For debugging, run flask app:

FLASK_APP=lib/postgap/Server.py FLASK_DEBUG=1 flask run

This by default serves on http://127.0.0.1:5000/ and gives an interactive debugger in your browser.
The debugger requires the pin code given in the log on the command line.
