(base) michael@ohm:~/gitrepos$ git clone https://github.com/lammps/lammps.git
Cloning into 'lammps'...
remote: Enumerating objects: 425797, done.
remote: Counting objects: 100% (259/259), done.
remote: Compressing objects: 100% (88/88), done.
remote: Total 425797 (delta 231), reused 176 (delta 171), pack-reused 425538 (from 3)
Receiving objects: 100% (425797/425797), 785.47 MiB | 6.67 MiB/s, done.
Resolving deltas: 100% (351163/351163), done.
Updating files: 100% (13833/13833), done.
(base) michael@ohm:~/gitrepos$ cd lammps/
(base) michael@ohm:~/gitrepos/lammps$ git fetch origin --tags
(base) michael@ohm:~/gitrepos/lammps$ git tag -l 'stable_2*' | sort -V | tail -n 5
stable_29Oct2020
stable_29Sep2021
stable_29Sep2021_update1
stable_29Sep2021_update2
stable_29Sep2021_update3