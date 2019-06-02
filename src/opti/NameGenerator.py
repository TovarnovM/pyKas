# %load OptiSwarm\NameGenerator.py
import random
import os

# from http://www.geocities.com/anvrill/names/cc_goth.html
PLACES = ['Adara', 'Adena', 'Adrianne', 'Alarice', 'Alvita', 'Amara', 'Ambika', 'Antonia', 'Araceli', 'Balandria',
          'Basha',
          'Beryl', 'Bryn', 'Callia', 'Caryssa', 'Cassandra', 'Casondrah', 'Chatha', 'Ciara', 'Cynara', 'Cytheria',
          'Dabria', 'Darcei',
          'Deandra', 'Deirdre', 'Delores', 'Desdomna', 'Devi', 'Dominique', 'Drucilla', 'Duvessa', 'Ebony', 'Fantine',
          'Fuscienne',
          'Gabi', 'Gallia', 'Hanna', 'Hedda', 'Jerica', 'Jetta', 'Joby', 'Kacila', 'Kagami', 'Kala', 'Kallie', 'Keelia',
          'Kerry',
          'Kerry-Ann', 'Kimberly', 'Killian', 'Kory', 'Lilith', 'Lucretia', 'Lysha', 'Mercedes', 'Mia', 'Maura',
          'Perdita', 'Quella',
          'Riona', 'Safiya', 'Salina', 'Severin', 'Sidonia', 'Sirena', 'Solita', 'Tempest', 'Thea', 'Treva', 'Trista',
          'Vala', 'Winta']


###############################################################################
# Markov Name model
# A random name generator, by Peter Corbett
# http://www.pick.ucam.org/~ptc24/mchain.html
# This script is hereby entered into the public domain
###############################################################################
class Mdict:
    def __init__(self):
        self.d = {}

    def __getitem__(self, key):
        if key in self.d:
            return self.d[key]
        else:
            raise KeyError(key)

    def add_key(self, prefix, suffix):
        if prefix in self.d:
            self.d[prefix].append(suffix)
        else:
            self.d[prefix] = [suffix]

    def get_suffix(self, prefix):
        l = self[prefix]
        return random.choice(l)


class NameGenerator:
    """
    A name from a Markov chain
    """

    def __init__(self, chainlen=2, namelst=PLACES, unique_garant=1000):
        """
        Building the dictionary
        """
        if chainlen > 10 or chainlen < 1:
            print("Chain length must be between 1 and 10, inclusive")
            chainlen = 5

        self.mcd = Mdict()
        oldnames = []
        self.chainlen = chainlen
        self.unique_garant = unique_garant
        self.namesets = set()

        for l in namelst:
            l = l.strip()
            oldnames.append(l)
            s = " " * chainlen + l
            for n in range(0, len(l)):
                self.mcd.add_key(s[n:n + chainlen], s[n + chainlen])
            self.mcd.add_key(s[len(l):len(l) + chainlen], "\n")

    def new(self, unique_set=None):
        """
        New name from the Markov chain
        """
        if unique_set:
            unique_set = set(unique_set)
        prefix = " " * self.chainlen
        name = ""
        suffix = ""
        while True:
            suffix = self.mcd.get_suffix(prefix)
            if suffix == "\n" or len(name) > 9:
                name = name.capitalize()
                if name in self.namesets:
                    prefix = " " * self.chainlen
                    name = ""
                    suffix = ""
                    continue
                if unique_set:
                    if name in unique_set:
                        prefix = " " * self.chainlen
                        name = ""
                        suffix = ""
                        continue
                if len(self.namesets) > self.unique_garant:
                    self.namesets = set()
                self.namesets.update([name])
                break
            else:
                name = name + suffix
                prefix = prefix[1:] + suffix
        return name

    @staticmethod
    def get_standart(chainlen=2, u_un=1e5):
        fname = os.path.dirname(__file__)+'/Names.txt'
        with open(fname, 'r', encoding='utf-8') as f:
            res = [l.strip() for l in f.readlines()]
        return NameGenerator(chainlen, res, u_un)


def main():
    from tqdm import tqdm
    fname = 'Names.txt'
    with open(fname, 'r', encoding='utf-8') as f:
        print(f'Открытие файла с именами {fname}')
        res = [l.strip() for l in f.readlines()]
    u_un = 100000
    print(f"Создание генератора имен с гарантированной генерацией {u_un} имен")
    ng = NameGenerator(2, res, u_un)

    n_count = u_un
    print(f"Создание списка из {n_count} имен")
    names = [ng.new() for i in tqdm(range(n_count))]
    uniquen = set(names)
    unique_raito = len(uniquen) / len(names)
    print(f'Процент уникальных имен среди {n_count} имен = {unique_raito*100}%')

    n_count = u_un * 10
    print(f"Создание списка из {n_count} имен")
    names = [ng.new() for i in tqdm(range(n_count))]
    uniquen = set(names)
    unique_raito = len(uniquen) / len(names)
    print(f'Процент уникальных имен среди {n_count} имен = {unique_raito*100}%')

    n_count = u_un * 100
    print(f"Создание списка из {n_count} имен")
    names = [ng.new() for i in tqdm(range(n_count))]
    uniquen = set(names)
    unique_raito = len(uniquen) / len(names)
    print(f'Процент уникальных имен среди {n_count} имен = {unique_raito*100}%')


if __name__ == '__main__':
    main()
