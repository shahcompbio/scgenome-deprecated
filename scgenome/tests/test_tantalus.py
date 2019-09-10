from scgenome import tantalus
from unittest import TestCase, main
from scgenome.TNode import TNode


class TestTantalus(TestCase):

    def test_spike_in(self):
        total_ncells = 100
        hmmcopy_tickets = ['SC-1936', 'SC-1937']
        sample_ids = [['SA921'], ['SA1090']]
        #hmmcopy_tickets = ['SC-1936']
        #sample_ids = [['SA921']]
        #hmmcopy_tickets = ['SC-1937']
        #sample_ids = [['SA1090']]

        results = tantalus.spike_in(total_ncells, hmmcopy_tickets, sample_ids,
                                    cached=True)

        print(results)