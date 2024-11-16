import argparse
import logging
from cadofactor import wudb
from . import ServerLauncher

serveraddress = '0.0.0.0'
serverport = 8080
db = None
threaded = True

parser = argparse.ArgumentParser('example server')
parser.add_argument("-d", "--debug", action='store_true',
                    default=False, help="enable debug messages")
parser.add_argument("-D", "--database", help='database uri')
parser.add_argument("-l", "--listen", default='0.0.0.0',
                    help='server address')
parser.add_argument("-p", "--port", default=8080, type=int,
                    help='server port')
parser.add_argument("-t", "--threads", default=None, type=int,
                    help='number of workers')
parser.add_argument("-c", "--cafile", default=None,
                    help='ssl certificate name')

args = parser.parse_args()

logger = logging.getLogger(__name__)
logger.setLevel(logging.NOTSET)

server = ServerLauncher(args.listen, args.port,
                        wudb.DBFactory(args.database),
                        threaded=args.threads,
                        debug=args.debug,
                        cafile=args.cafile,
                        uploaddir='/tmp/uploads',
                        whitelist=['0.0.0.0/0'],
                        nrsubdir=4)
server.serve()
