from Bio import SeqIO
from random import (random, randint, choice)
from enum import Enum, auto


class Event(Enum):
    INSERTION = auto()
    DELETION = auto()
    SUBSTITUTION = auto()
    STAY = auto()


def rnd_nuc() -> str:
    return choice("ATGC")


def rnd_read(ref_seq: str, read_len: int) -> zip:
    pos = randint(0, len(ref_seq) - 1 - read_len)
    return zip(range(pos, pos + read_len), ref[pos: pos + read_len])


def make_event(event: Event, nuc: str) -> str:
    if event is Event.INSERTION:
        return nuc + rnd_nuc()
    if event is Event.DELETION:
        return ""
    if event is Event.SUBSTITUTION:
        return rnd_nuc()
    if event is Event.STAY:
        return nuc


if __name__ == "__main__":

    readAmount = 10 ** 4
    readLen = 200
    allEvents = list(Event.__members__.values())
    ref = next(SeqIO.parse("./ref.fa", "fasta")).seq

    # reads only with substitutions
    sub_prob = 0.01
    events = {
        p: Event.SUBSTITUTION if random() < sub_prob else Event.STAY
        for p in range(len(ref))
    }
    with open("./reads/reads_substitution.fasta", "w") as reads:
        for idx in range(readAmount):
            seq = rnd_read(ref, readLen)
            readSeq = "".join(
                [make_event(events[pos], nuc) for pos, nuc in seq]
            )
            reads.write(">{}\n{}\n".format(idx, readSeq))

    # reads only with any types of event
    event_prob = 0.05
    events = {
        p: choice([Event.SUBSTITUTION, Event.DELETION, Event.INSERTION])
        if random() < sub_prob else Event.STAY
        for p in range(len(ref))
    }
    with open("./reads/reads_all_events.fasta", "w") as reads:
        for idx in range(readAmount):
            seq = rnd_read(ref, readLen)
            readSeq = "".join(
                [make_event(events[pos], nuc) for pos, nuc in seq]
            )
            readSeq = readSeq[:readLen]
            reads.write(">{}\n{}\n".format(idx, readSeq))
