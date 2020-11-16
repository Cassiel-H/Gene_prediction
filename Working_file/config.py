# Jiachen Huo
# 260779330


from collections import defaultdict
from Q1a import Q1a


class config:
    def __init__(self):
        a = Q1a()
        self.initial = [1, 0, 0, 0]
        self.transition = Transition(a)
        self.emission = Emission(a)
        self.state = {'Inter': 1, 'Start': 3, 'Middle': '3', 'Stop': 3}


def Emission(a):
    inter = a.nuc_frequency
    start = a.startTable.set_index('codon')['frequency'].to_dict()
    middle = a.middleCodon.set_index('codon')['frequency'].to_dict()
    stop = a.stopTable.set_index('codon')['frequency'].to_dict()
    e_table = [
        defaultdict(float, inter),
        defaultdict(float, start),
        defaultdict(float, middle),
        defaultdict(float, stop)
    ]
    return e_table


def Transition(a):
    interLength = a.avgInterLen
    genLength = a.avgenicLen / 3
    IM = IStop = StartI = StartStart = StartStop = MI = MStart = StopStart = StopM = StopStop = 0
    StartM = StopI = 1
    II = ((interLength - 1) / interLength)
    IStart = (1 / interLength)
    MM = ((genLength - 1) / genLength)
    MStop = (1 / genLength)
    t_table = [
        [II, IStart, IM, IStop],
        [StartI, StartStart, StartM, StartStop],
        [MI, MStart, MM, MStop],
        [StopI, StopStart, StopM, StopStop]
    ]
    # print(t_table)
    return t_table


def main():
    fig = config()


if __name__ == '__main__':
    main()
