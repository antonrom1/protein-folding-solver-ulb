# cours informatique fondamentale 2021-2022
# PROJET: repliage de proteines

# necessite l'installation de la librairie PySAT et de la librairie func_timeout
import itertools
import logging
import multiprocessing
import sys
import time
from optparse import OptionParser
from pprint import pp

from pysat.card import *

from pysat.solvers import Glucose4

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

# on doit TOUJOURS donner une sequence
# * lorsqu'une borne est donnee, votre programme doit tester que le meilleur score de la sequence est superieur ou egale a cette borne
# * lorsqu'aucune borne n'est donnee, alors votre programme doit calculer le meilleur score pour la sequence, par defaut en utilisant une recherche par dichotomie, et en utilisant une methode incrementale si l'option -i est active
#
# l'option -v vous permet de creer un mode 'verbose'
# si l'option -t est active, alors le code execute uniquement la fonction test_code() implementee ci-dessous, qui vous permet de tester votre code avec des exemples deja fournis. Si l'execution d'un test prend plus que TIMEOUT secondes (fixe a 10s ci-dessous), alors le test s'arrete et la fonction passe au test suivant

parser = OptionParser()
parser.add_option("-s", "--sequence", dest="seq", action="store",
                  help="specify the input sequence")
parser.add_option("-b", "--bound", dest="bound", action="store",
                  help="specify a lower bound on the score", type="int")
parser.add_option("-p", "--print", dest="affichage_sol", action="store_true",
                  help="print solution", default=False)
parser.add_option("-i", "--incremental", dest="incremental", action="store_true",
                  help="incremental mode: try small bounds first and increment", default=False)
parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                  help="verbose mode", default=False)
parser.add_option("-t", "--test", dest="test", action="store_true",
                  help="testing mode", default=False)
parser.add_option("-d", "--debug", dest="debug", action="store_true",
                  help="debug logging", default=False)

(options, args) = parser.parse_args()

affichage_sol = options.affichage_sol
verb = options.verbose
incremental = options.incremental
test = options.test
idp = IDPool()
logging.basicConfig(level=logging.DEBUG if options.debug else (logging.INFO if options.verbose else logging.WARNING))

###############################################################################################


def max_bound_for_seq_of_length(n):
    from math import sqrt
    from math import ceil
    return ceil(2 * (sqrt(n) - 1) * sqrt(n))




def exist_sol(seq, bound):
    if not seq.isnumeric():
        raise ValueError("Sequence must be numeric")

    # eliminate simple cases
    if bound > max_bound_for_seq_of_length(len(seq)):
        return False

    # trivial case of only '1' in sequence
    if all(x == '1' for x in seq) and bound < max_bound_for_seq_of_length(len(seq)):
        return True

    cnf = CNF()
    indices = set(range(len(seq)))
    possible_matrix_values = indices | {None, }


    logging.info("check if every matrix position has only one value")
    for i, j in itertools.product(range(len(seq)), repeat=2):
        cnf.extend(CardEnc.equals([idp.id((i, j, v)) for v in possible_matrix_values], 1, vpool=idp))

    only_one_value_per_position_clauses_count = len(cnf.clauses)
    logging.debug(f'Generated {only_one_value_per_position_clauses_count} clauses for "only one value per position"')


    logging.info("check if every value (except None) appears only once")
    for v in possible_matrix_values - {None}:
        cnf.extend(CardEnc.equals([idp.id((i, j, v)) for i, j in itertools.product(range(len(seq)), repeat=2)], 1, vpool=idp))

    only_one_index_in_matrix_clauses_count = len(cnf.clauses) - only_one_value_per_position_clauses_count
    logging.debug(f'Generated {only_one_index_in_matrix_clauses_count} clauses for "only one index in matrix"')


    logging.info("check if every value is adjacent to its neighbours")
    # check whether i is adjacent to i - 1 and i + 1
    for i, j in itertools.product(range(len(seq)), repeat=2):
        # TODO: this can be optimized by skipping 1/2 of the values because they are already checked
        for v in range(1, len(seq) - 1):
            # [(i, j, v) => [(i - 1, j, v - 1) | (i + 1, j, v - 1) | (i, j - 1, v - 1) | (i, j + 1, v - 1)]]
            # a => (b | c | d | e) <=> !a | b | c | d | e

            v1p1_clause = [-idp.id((i, j, v))]
            v1m1_clause = [-idp.id((i, j, v))]

            for i2, j2 in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                if 0 <= i2 < len(seq) and 0 <= j2 < len(seq):
                    v1p1_clause.append(idp.id((i2, j2, v - 1)))
                    v1m1_clause.append(idp.id((i2, j2, v + 1)))
            cnf.append(v1m1_clause)
            cnf.append(v1p1_clause)

    adjacent_values_clauses_count = len(cnf.clauses) - only_one_value_per_position_clauses_count - only_one_index_in_matrix_clauses_count
    logging.debug(f'Generated {adjacent_values_clauses_count} clauses for "adjacent values"')

    logging.info("check if there are at least 'bound' adjacent 1s")
    # we have to user card.atleast here
    clauses_for_atleast = []
    for i, j in itertools.product(range(len(seq)), repeat=2):
        # indices of the 1s in the sequence
        if (i + j) % 2 == 0:
            continue
        indices_of_1s = [i for i, v in enumerate(seq) if v == '1']
        for v1, v2 in itertools.combinations(indices_of_1s, 2):
            for i2, j2 in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                if not (0 <= i2 < len(seq) and 0 <= j2 < len(seq)):
                    continue
                a = idp.id((i, j, v1))
                b = idp.id((i2, j2, v2))
                c = idp.id(((i, i2), (j, j2), (v1, v2)))
                # append (a & b) <=> c
                # (¬c ∨ a) ∧ (¬c ∨ b) ∧ (¬a ∨ ¬b ∨ c)
                cnf.extend([[-c, a], [-c, b], [a, b, -c]])
                clauses_for_atleast.append(c)

    if bound > len(clauses_for_atleast):
        return False

    cnf.extend(CardEnc.atleast(clauses_for_atleast, bound, vpool=idp))

    atleast_clauses_count = len(cnf.clauses) - only_one_value_per_position_clauses_count - only_one_index_in_matrix_clauses_count - adjacent_values_clauses_count
    logging.debug(f'Generated {atleast_clauses_count} clauses for "atleast"')

    logging.info(f"A total of {len(cnf.clauses)} clauses were generated")

    # we solve the problem
    logging.debug("solving...")
    s = Glucose4(use_timer=True)
    s.append_formula(cnf.clauses, no_return=True)
    s.solve()
    logging.debug("done")

    return s.get_model() is not None


def compute_max_score(
        seq):  # calcul le meilleur score pour la sequence seq, il doit donc retourne un entier, methode utilisee: dichotomie par defaut, incrementale si l'option -i est active
    # A COMPLETER
    return (10)  # a modifier


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
    # maxlength = 25
    # m = 2
    # l = []
    # for n in range(maxlength):
    # for i in range(m):
    # seq = gen(n+2)
    # print(seq)
    # b = compute_max_score(seq)
    # print(int(b))
    # l.append((seq,b))
    # print(l)

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
        print("DEBUT DU TEST DE SATISFIABILITE")
        if exist_sol(options.seq, options.bound):
            print("SAT")
            # print(model)
            # model_set = set(model)
            # output_table = [['-'] * len(seq) for _ in range(len(seq))]
            # for i, j in itertools.product(range(len(seq)), repeat=2):
            #     for v in range(len(seq)):
            #         if idp.id((i, j, v)) in model_set:
            #             if output_table[i][j] != '-':
            #                 print("ERROR: multiple values for cell (" + str(i) + "," + str(j) + ")")
            #                 exit(1)
            #             output_table[i][j] = seq[v]
            # print(*output_table, sep='\n')
            #
            # cs = []
            # for i, j in itertools.product(range(len(seq)), repeat=2):
            #     # indices of the 1s in the sequence
            #     indices_of_1s = [i for i, v in enumerate(seq) if v == '1']
            #     for v1, v2 in itertools.combinations(indices_of_1s, 2):
            #         for i2, j2 in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
            #             if not (0 <= i2 < len(seq) and 0 <= j2 < len(seq)):
            #                 continue
            #             a = idp.id((i, j, v1))
            #             b = idp.id((i + i2, j + j2, v2))
            #             i1i2_sorted = (i, i2)# if i < i2 else (i2, i)
            #             j1j2_sorted = (j, j2)# if j < j2 else (j2, j)
            #             v1v2_sorted = (v1, v2)# if v1 < v2 else (v2, v1)
            #             c = idp.id((i1i2_sorted, j1j2_sorted, v1v2_sorted))
            #             cs.append(c)
            # print(len(cs), "cs in total")
            #
            # num_cs_true = 0
            # for c in cs:
            #     if c in model_set:
            #         num_cs_true += 1
            #
            # print(num_cs_true, "cs are true")



        else:
            print("UNSAT")
        print("FIN DU TEST DE SATISFIABILITE")

    elif not (incremental):
        # on affiche le score maximal qu'on calcule par dichotomie
        # on affiche egalement un plongement de score maximal si l'option d'affichage est active
        print("DEBUT DU CALCUL DU MEILLEUR SCORE PAR DICHOTOMIE")
        # A COMPLETER
        print("FIN DU CALCUL DU MEILLEUR SCORE")

    elif not test:
        # Pareil que dans le cas precedent mais avec la methode incrementale
        # A COMPLETER
        print("DEBUT DU CALCUL DU MEILLEUR SCORE PAR METHODE INCREMENTALE")
        # A COMPLETER
        print("FIN DU CALCUL DU MEILLEUR SCORE")
