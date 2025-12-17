# commit_for_second_exercise
Τα test μου στην αναφορά έχουν γίνει στο rome partition 
Οδηγίες για να τρέξει το πρόγραμμα:
cd $SCRATCH Μπαίνουμε στον SCRATCH folder  

wget https://suitesparse-collection-website.herokuapp.com/mat/SNAP/com-Friendster.mat   
π.χ για τον com-Friendster δειχνω εδω . Κανουμε wget το .mat link απο το suitesparse για τον γραφο που θελουμε

cd $HOME ή cd $HOME/cloned_repository_folder   μπαινουμε πισω εκει που εχουμε κανει clone τα αρχεια 

sbatch submission-script.sh com-Friendster.mat       τρεχουμε το script με ορισμα το .mat αρχειο του γραφου που θελουμε να τρεξουμε

Την πρωτη φορα που κανουμε τεστ γενικα, σεταρεται στο SCRATCH ενα environment ωστε να τρεξει το python script μου.
Την πρωτη φορα που κανουμε test για ενα συγκεκριμενο γραφο τωρα, τρεχει το script και βγαζει τον φακελο με τα binary αρχει και το .txt
Οποτε εχουμε μια παραπανω καθυστερση αυτες τις φορες. 1 στο πρωτο τεστ απο ολα και 1 για καθε διαφορετικο γραφο
Οσο εχει να κανει με τις παραμετρους του slurm τις οριζουμε αλλαζοντας το script οπως επιθυμουμε
Στα δικα μου τεστ ετρεξα με 100G τον com_Friendster και με 60G τον Kmer_V1r
