from __future__ import absolute_import, division
from dlstbx.util.colorstreamhandler import ColorStreamHandler
import json
import logging
from optparse import OptionParser, SUPPRESS_HELP
import sys
import time
from workflows.transport.stomp_transport import StompTransport



class Worker(object):

  def __init__(self, wait_time=5):

    self.wait_time = wait_time
    self.heartbeat_timestamp = None

    self.transport = StompTransport()
    self.transport.connect()
    self.transport.subscribe_broadcast("heartbeat", self.read_heartbeat)
    self.transport.subscribe("outbound", self.read)

    self.wait()

  def timeout(self):
    from time import time
    if self.heartbeat_timestamp is None:
      self.heartbeat_timestamp = time()
    else:
      if time() - self.heartbeat_timestamp > self.wait_time:
        return True
    return False

  def wait(self):
    from time import time, sleep
    try:
      while not self.timeout():
        sleep(1)
    except KeyboardInterrupt:
      pass

  def send(self, message):
    print "Sending"
    self.transport.send(
      'inbound',
      message
    )

  def read_heartbeat(self, header, message):
    from time import time
    self.heartbeat_timestamp = time()

  def read(self, header, message):
      message = {
        "result" : "Hey Hey!"
      }
      self.send(message)



if __name__ == '__main__':
  configuration = '/dls_sw/apps/zocalo/secrets/credentials-james.cfg'

  from workflows.transport.stomp_transport import StompTransport

  StompTransport.load_configuration_file(configuration)

  worker = Worker()
