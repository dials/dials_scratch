from __future__ import division
import workflows
from workflows.transport.stomp_transport import StompTransport
from threading import Thread, Event


class Sender(object):

  def __init__(self, transport):
    self.transport = transport

  def send(self, message):
    self.transport.send(
      'outbound',
      message
    )

class Receiver(object):

  def __init__(self, transport):
    self.transport = transport

    self.message_list = []

    self.transport.subscribe("inbound", self.read)

  def read(self, header, message):
    self.message_list.append(message)


class HeartBeat(Thread):

  def __init__(self, transport):
    super(HeartBeat, self).__init__()
    self.transport = transport
    self.stop_event = Event()

  def run(self):
    while not self.stop_event.wait(0.5):
      self.transport.broadcast("heartbeat", "")

  def stop(self):
    self.stop_event.set()


class Master(object):

  def __init__(self):

    self.transport = StompTransport()
    self.transport.connect()

    self.heartbeat = HeartBeat(self.transport)
    self.heartbeat.start()

    self.sender = Sender(self.transport)
    self.receiver = Receiver(self.transport)
    self.num_messages = 0

  def send(self, message):
    self.num_messages += 1
    self.sender.send(message)

  def receive(self):

    try:
      from time import sleep
      while len(self.receiver.message_list) != self.num_messages:
        print len(self.receiver.message_list)
        sleep(1)
      self.heartbeat.stop()
    except KeyboardInterrupt:
      self.heartbeat.stop()
      raise

    return self.receiver.message_list


def cluster_map(
    func,
    iterable,
    callback=None,
    nslots=1,
    njobs=1,
    job_category="low"):
  from multiprocessing import Process

  processes = [Process(target=Worker()) for i in range(len(njobs))]
  for p in processes:
    p.start()




  for p in processes:
    p.join()


class Data(object):

  def __init__(self, a, b, c):
    self.a = a
    self.b = b
    self.c = c


if __name__ == '__main__':
  import cPickle as pickle

  configuration = '/dls_sw/apps/zocalo/secrets/credentials-james.cfg'

  StompTransport.load_configuration_file(configuration)

  data = Data(1, 2, 3)

  message = {
    'hello': 'world',
    'data' : pickle.dumps(data),
  }

  master = Master()
  master.send(message)
  master.send(message)
  master.send(message)
  master.send(message)
  print master.receive()


  print "\nSubmitted."
