import threading


class Pauser(object):
    def __init__(self):
        self.cv = threading.Condition()
        self.paused = False

    def pause(self):
        self.paused = True

    def play(self):
        self.paused = False
        with self.cv:
            self.cv.notify_all()

    def kpp(self):
        if self.paused:
            with self.cv:
                self.cv.wait()


if __name__ == '__main__':
    import time
    import logging

    logging.basicConfig(level=logging.DEBUG,
                        format='(%(threadName)-9s) %(message)s', )


    def consumer(pser, doner):
        logging.debug('Consumer thread started ...')
        while not doner.paused:
            logging.debug('Consumer thread working ...')
            time.sleep(2)
            pser.kpp()


    pser = Pauser()
    doner = Pauser()
    cs1 = threading.Thread(name='consumer1', target=consumer, args=(pser, doner))
    cs2 = threading.Thread(name='consumer2', target=consumer, args=(pser, doner))

    cs1.start()
    time.sleep(1)
    cs2.start()
    time.sleep(5)
    pser.pause()
    logging.debug('Pause ...')
    time.sleep(5)
    pser.play()
    logging.debug('Playing ...')
    time.sleep(5)
    doner.pause()
