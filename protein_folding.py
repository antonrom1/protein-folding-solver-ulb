# cours informatique fondamentale 2021-2022
# PROJET: repliage de proteines

# necessite l'installation de la librairie PySAT et de la librairie func_timeout
import logging
import multiprocessing
import sys
import time
from optparse import OptionParser

from pysat.card import *

# on doit TOUJOURS donner une sequence
# * lorsqu'une borne est donnee, votre programme doit tester que le meilleur score de la sequence est superieur ou egale a cette borne
# * lorsqu'aucune borne n'est donnee, alors votre programme doit calculer le meilleur score pour la sequence, par defaut en utilisant une recherche par dichotomie, et en utilisant une methode incrementale si l'option -i est active
#
# l'option -v vous permet de creer un mode 'verbose'
# si l'option -t est active, alors le code execute uniquement la fonction test_code() implementee ci-dessous, qui vous permet de tester votre code avec des exemples deja fournis. Si l'execution d'un test prend plus que TIMEOUT secondes (fixe a 10s ci-dessous), alors le test s'arrete et la fonction passe au test suivant
from src.folding import FoldingSolver, FoldingSolution, FoldedProtein
from src.sequence import Sequence
from src.utils import get_max_score_bound_for_length

##### OPTIONS POUR L'UTILISATION EN LIGNE DE COMMANDE ###############
# Usage: folder.py [options]
# Options:
# -h, --help            show this help message and exit
# -s SEQ, --sequence=SEQ
# specify the input sequence
# -b BOUND, --bound=BOUND
# specify a lower bound on the score
# -p, --print           print solution
# -i, --incremental     incremental mode: try small bounds first and increment
# -v, --verbose         verbose mode
# -t, --test            testing mode

parser = OptionParser()
parser.add_option("-s", "--sequence", dest="seq", action="store",
                  help="specify the input sequence")
parser.add_option("-b", "--bound", dest="bound", action="store",
                  help="specify a lower bound on the score", type="int")
parser.add_option("-p", "--print", dest="affichage_sol", action="store_true",
                  help="print solution", default=False)
parser.add_option("-i", "--incremental", dest="incremental", action="store_true",
                  help="incremental mode: try small bounds first and increment", default=False)
parser.add_option("-v", "--verbose", dest="loglevel", action="store_const",
                  help="verbose mode", const=logging.INFO, default=logging.WARNING)
parser.add_option("-t", "--test", dest="test", action="store_true",
                  help="testing mode", default=False)
parser.add_option("-d", "--debug", dest="loglevel", action="store_const",
                  help="debug mode", const=logging.DEBUG)

(options, args) = parser.parse_args()

affichage_sol = options.affichage_sol
incremental = options.incremental
test = options.test
idp = IDPool()
logging.basicConfig(level=options.loglevel)
logger = logging.getLogger(__name__)


###############################################################################################


def exist_sol(seq, bound):
    """return True if there is a solution with score >= bound"""
    sequence = Sequence(seq)
    sol = FoldingSolver(sequence).solve(bound)
    if options.affichage_sol:
        if sol.is_sat:
            print(sol.solution)
        else:
            print("No solution found")
    return sol.is_sat


def compute_max_score(seq):
    seq = Sequence(seq)
    solver = FoldingSolver(seq)
    min_score = solver.sequence.get_flat_sequence_score()
    if incremental:
        max_bound = get_max_score_bound_for_length(seq.number_of_1s)
        solutions = []

        for i in range(min_score, max_bound + 1):
            sol = solver.solve(i)
            if sol.is_sat:
                solutions.append(sol)
            else:
                break

        solution = solutions[-1]

        if options.affichage_sol:
            print(solution.solution)

        return solution.bound

    else:
        max_bound = get_max_score_bound_for_length(seq.number_of_1s)
        min_bound = min_score
        solution = FoldingSolution(sequence=seq, bound=0, score=min_bound,
                                   solution=FoldedProtein.from_straight_sequence(seq))
        while min_bound < max_bound:
            bound = (min_bound + max_bound) // 2
            res = solver.solve(bound)
            if res.is_sat:
                solution = res
                min_bound = bound + 1
            else:
                max_bound = bound

        previous_solution = solution

        solution = solver.solve(min_bound)

        if solution.is_unsat:
            solution = previous_solution

        if options.affichage_sol:
            if solution.is_sat:
                print(solution.solution)
            else:
                logger.error("No solution found")
                return 0

        return solution.bound


####################################################################
########### CE CODE NE DOIT PAS ETRE MODIFIE #######################
####################################################################
def worker_exist_sol(queue, seq, bound):
    """worker function"""
    ret = queue.get()
    try:
        s = exist_sol(seq, bound)
        ret["exist_sol"] = s
    except Exception as e:
        ret["exist_sol"] = e
    queue.put(ret)


def worker_max_score(queue, seq):
    """worker function"""
    ret = queue.get()
    try:
        s = compute_max_score(seq)
        ret["exist_sol"] = s
    except Exception as e:
        ret["exist_sol"] = e
    queue.put(ret)


def test_code():
    examples = [('00', 0), ('1', 0), ('01000', 0), ('00110000', 1), ('11', 1), ('111', 2), ('1111', 4), ('1111111', 8),
                ("111111111111111", 22), ("1011011011", 7), ("011010111110011", 12), ("01101011111000101", 11),
                ("0110111001000101", 8), ("000000000111000000110000000", 5), ('100010100', 0),
                ('01101011111110111', 17), ('10', 0), ('10', 0), ('001', 0), ('000', 0), ('1001', 1), ('1111', 4),
                ('00111', 2), ('01001', 1), ('111010', 3), ('110110', 3), ('0010110', 2), ('0000001', 0),
                ('01101000', 2), ('10011111', 7), ('011001101', 5), ('000110111', 5), ('0011000010', 2),
                ('1000010100', 2), ('11000111000', 5), ('01010101110', 4), ('011001100010', 5), ('010011100010', 5),
                ('1110000110011', 8), ('1000101110001', 4), ('11010101011110', 10), ('01000101000101', 0),
                ('111011100100000', 8), ('000001100111010', 6), ('0110111110011000', 11), ('0011110010110110', 11),
                ('01111100010010101', 11), ('10011011011100101', 12), ('101111101100101001', 13),
                ('110101011010101010', 9), ('1111101010000111001', 14), ('0111000101001000111', 11),
                ('10111110100001010010', 12), ('10110011010010001110', 11)]

    TIMEOUT = 10
    ret = {"ret_val": None}  # to store the returned value of process calls

    # SAT TESTS
    total_sat_tests = 0
    sat_tests_success = 0
    sat_tests_failure = 0
    timeouts_sat_tests = 0
    exceptions_sat_tests = 0

    # UNSAT TESTS
    total_unsat_tests = 0
    total_unsat_test = 0
    unsat_tests_success = 0
    unsat_tests_failure = 0
    timeouts_unsat_tests = 0
    exceptions_unsat_tests = 0

    # MAXSCORES TEST
    correct_maxscores = 0
    incorrect_maxscores = 0
    total_maxscores = 0
    timeouts_maxscores = 0
    exceptions_maxscores = 0

    # sur cet ensemble de tests, votre methode devrait toujours retourner qu'il existe une solution
    print("\n****** Test de satisfiabilite ******\n")
    for (seq, maxbound) in examples:
        # initialize queue
        queue = multiprocessing.Queue()
        queue.put(ret)
        ret = queue.get()
        ret["exist_sol"] = None
        queue.put(ret)

        # increase counter
        total_sat_tests += 1
        # set bound
        bound = int(maxbound / 2)
        print("sequence: " + seq + " borne: " + str(bound), end='')
        sys.stdout.flush()

        p = multiprocessing.Process(target=worker_exist_sol, args=(queue, seq, bound))
        p.start()
        p.join(timeout=TIMEOUT)

        if p.exitcode is None:
            print(" ---> timeout")
            timeouts_sat_tests += 1
            p.terminate()
        else:
            r = queue.get()["exist_sol"]
            if isinstance(r, Exception):
                exceptions_sat_tests += 1
                print(" ---> exception raised")
            else:
                if r == True:
                    sat_tests_success += 1
                    print(" ---> succes")
                else:
                    sat_tests_failure += 1
                    print(" ---> echec")

    # sur cet ensemble de tests, votre methode devrait toujours retourner qu'il n'existe pas de solution
    print("\n****** Test de d'insatisfiabilite ******\n")
    for (seq, maxbound) in examples:
        # initialize queue
        queue = multiprocessing.Queue()
        queue.put(ret)
        ret = queue.get()
        ret["exist_sol"] = None
        queue.put(ret)

        # increase counter
        total_unsat_tests += 1
        # set bound
        bound = maxbound + 1
        print("sequence: " + seq + " borne: " + str(bound), end='')
        sys.stdout.flush()

        p = multiprocessing.Process(target=worker_exist_sol, args=(queue, seq, bound))
        p.start()
        p.join(timeout=TIMEOUT)

        if p.exitcode is None:
            print(" ---> timeout")
            timeouts_unsat_tests += 1
            p.terminate()
        else:
            r = queue.get()["exist_sol"]
            if isinstance(r, Exception):
                exceptions_unsat_tests += 1
                print(" ---> exception raised")
            else:
                if r == True:
                    unsat_tests_failure += 1
                    print(" ---> echec")
                else:
                    unsat_tests_success += 1
                    print(" ---> succes")

    # sur cet ensemble de tests, votre methode devrait retourner le meilleur score. Vous pouvez utiliser la methode par dichotomie ou incrementale, au choix
    print("\n****** Test de calcul du meilleur score ******\n")
    for (seq, maxbound) in examples:
        # initialize queue
        queue = multiprocessing.Queue()
        queue.put(ret)
        ret = queue.get()
        ret["exist_sol"] = None
        queue.put(ret)

        total_maxscores += 1
        print("sequence: " + seq + " borne attendue: " + str(maxbound), end='')
        sys.stdout.flush()

        p = multiprocessing.Process(target=worker_max_score, args=(queue, seq,))
        p.start()
        p.join(timeout=TIMEOUT)

        if p.exitcode is None:
            print(" ---> timeout")
            timeouts_maxscores += 1
            p.terminate()
        else:
            r = queue.get()["exist_sol"]
            if isinstance(r, Exception):
                exceptions_maxscores += 1
                print(r)
                print(" ---> exception raised")
            else:
                if r != maxbound:
                    incorrect_maxscores += 1
                    print(" ---> echec (borne retournee est " + str(r) + ")")
                else:
                    correct_maxscores += 1
                    print(" ---> succes")

    print("\nRESULTATS TESTS\n")

    print("Nombre de total de tests de satisfiabilite : " + str(total_sat_tests))
    print("      Correctement repondues: " + str(sat_tests_success))
    print("      Incorrectement repondues: " + str(sat_tests_failure))
    print("      Timeouts: " + str(timeouts_sat_tests))
    print("      Exceptions: " + str(exceptions_sat_tests) + "\n")

    print("Nombre de total de tests d' insatisfiabilite : " + str(total_unsat_tests))
    print("      Correctement repondues: " + str(unsat_tests_success))
    print("      Incorrectement repondues: " + str(unsat_tests_failure))
    print("      Timeouts: " + str(timeouts_unsat_tests))
    print("      Exceptions: " + str(exceptions_unsat_tests) + "\n")

    print("Nombre de total de tests de calcul du meilleur score : " + str(total_maxscores))
    print("      Correctement repondues: " + str(correct_maxscores))
    print("      Incorrectement repondues: " + str(incorrect_maxscores))
    print("      Timeouts: " + str(timeouts_maxscores))
    print("      Exceptions: " + str(exceptions_maxscores) + "\n")


##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################

if __name__ == '__main__':

    if test:
        print("Let's test your code")
        start_time_test = time.time()
        test_code()
        end_time_test = time.time()
        print("Temps total pour le test: ", end_time_test - start_time_test, "sec")

    elif options.bound != None:
        # cas ou la borne est fournie en entree: on test si la sequence (qui doit etre donnee en entree) a un score superieur ou egal a la borne donnee
        # si oui, on affiche "SAT". Si l'option d'affichage est active, alors il faut egalement afficher une solution
        seq = options.seq
        logger.debug(f"Searching for a solution for sequence {seq} with bound {options.bound}")
        if exist_sol(seq, options.bound):
            print("SAT")
        else:
            print("UNSAT")

    elif not (incremental):
        # on affiche le score maximal qu'on calcule par dichotomie
        # on affiche egalement un plongement de score maximal si l'option d'affichage est active
        print("DEBUT DU CALCUL DU MEILLEUR SCORE PAR DICHOTOMIE")
        start = time.perf_counter()
        print(compute_max_score(options.seq))
        end = time.perf_counter()

        print(f"Temps d'execution : {end - start}")
        print("FIN DU CALCUL DU MEILLEUR SCORE")

    elif not test:
        # Pareil que dans le cas precedent mais avec la methode incrementale
        # A COMPLETER
        print("DEBUT DU CALCUL DU MEILLEUR SCORE PAR METHODE INCREMENTALE")
        print(compute_max_score(options.seq))
        print("FIN DU CALCUL DU MEILLEUR SCORE")
