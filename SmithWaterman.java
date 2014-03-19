
/*
 * http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
 */

public class SmithWaterman {

    // >sp|P51690|ARSE_HUMAN Arylsulfatase E OS=Homo sapiens GN=ARSE PE=1 SV=2
    // http://www.uniprot.org/uniprot/P51690.fasta

    static String protein1 =
            "MLHLHHSCLCFRSWLPAMLAVLLSLAPSASSDISASRPNILLLMADDLGIGDIGCYGNNT" +
                    "MRTPNIDRLAEDGVKLTQHISAASLCTPSRAAFLTGRYPVRSGMVSSIGYRVLQWTGASG" +
                    "GLPTNETTFAKILKEKGYATGLIGKWHLGLNCESASDHCHHPLHHGFDHFYGMPFSLMGD" +
                    "CARWELSEKRVNLEQKLNFLFQVLALVALTLVAGKLTHLIPVSWMPVIWSALSAVLLLAS" +
                    "SYFVGALIVHADCFLMRNHTITEQPMCFQRTTPLILQEVASFLKRNKHGPFLLFVSFLHV" +
                    "HIPLITMENFLGKSLHGLYGDNVEEMDWMVGRILDTLDVEGLSNSTLIYFTSDHGGSLEN" +
                    "QLGNTQYGGWNGIYKGGKGMGGWEGGIRVPGIFRWPGVLPAGRVIGEPTSLMDVFPTVVR" +
                    "LAGGEVPQDRVIDGQDLLPLLLGTAQHSDHEFLMHYCERFLHAARWHQRDRGTMWKVHFV" +
                    "TPVFQPEGAGACYGRKVCPCFGEKVVHHDPPLLFDLSRDPSETHILTPASEPVFYQVMER" +
                    "VQQAVWEHQRTLSPVPLQLDRLGNIWRPWLQPCCGPFPLCWCLREDDPQ";
    // >tr|G3RJQ1|G3RJQ1_GORGO Uncharacterized protein OS=Gorilla gorilla gorilla GN=101154474 PE=4 SV=1
    // http://www.uniprot.org/uniprot/G3RJQ1.fasta
    static String protein2 =
            "MLHLHHSCLCFRSWLAAMLAVLLSLAPSASSDISASRPNILLLMADDLGIGDIGCYGNNT" +
                    "MRTPNIDRLAEDGVKLTQHISAASLCTPSRAAFLTGRYPVRSGMVSSIGYRVLQWTGASG" +
                    "GLPTNETTFAKILKEKGYATGLIGKWHLGLNCESASDHCHHPLHHGFDHFYGMPFSLMGD" +
                    "CARWELSEKRVNLEQKLNFLFQVLALVALTLVAGKLTHLIPVSWMPVIWSALSAVLLLAS" +
                    "SYFVGALIVHADCFLMRNHTITEQPMCFQRTTPLILQEVASFLKRNKHGPFLLFVSFLHV" +
                    "HIPLITMENFLGKSLHGLYGDNVEEMDWMVGRILDTLDVEGLSNSTLIYFTSDHGGSLEN" +
                   // "QLGNTQYGGWNGIYKDTGGKGMGGWEGGIRVPGIFRWPGVLPAGRVIGEPTSLMDVFPTV" +
                    "VRLAGGEVPQDRVIDGQDLLPLLLGTAQHSDHEFLMHYCERFLHAARWHQRDRGTMWKVH" +
                    "FVTPVFQPEGAGACYGRKVCPCFGEKVVHHDPPLLFDLSRDPSETHILTPASEPMFYQVM" +
                    "ERVQQAVREHQRTLSPVPLQLDRLGNIWRPWLQPCCGPFPLCWCLREDDPQ";

    public static void main(String[] args) {

        double score1 = new SmithWaterman().apply(protein1, protein1);
        double score2 = new SmithWaterman().apply(protein1, protein2);
        double score3 = new SmithWaterman().apply(protein2, protein2);

        System.out.println("score1 = " + score1);
        System.out.println("score2 = " + score2);
        System.out.println("score3 = " + score3);

    }

    private double gapOpeningPenalty;
    private double gapExtensionPenalty;
    public SmithWaterman() {
        this(8,0);
    }


    public SmithWaterman(double gapOpeningPenalty, double gapExtensionPenalty) {
        this.gapOpeningPenalty = gapOpeningPenalty;
        this.gapExtensionPenalty = gapExtensionPenalty;
    }

    public double apply(String s1, String s2) {

        double score = 0;
        char[] s1Chars = s1.toUpperCase().toCharArray();
        char[] s2Chars = s2.toUpperCase().toCharArray();

        int s1Length = s1.length() + 1;
        int s2Length = s2.length() + 1;

        double f; // score of alignment x1...xi to y1...yi if xi aligns to yi
        double[] s1Scores = new double[s2Length]; // score if xi aligns to a gap after yi
        double h; // score if yi aligns to a gap after xi
        double[] s2Scores = new double[s2Length]; // best score of alignment x1...xi to
        // y1...yi
        double vDiagonal;

        s1Scores[0] = Double.NEGATIVE_INFINITY;
        s2Scores[0] = 0;

        init(s2Length, s1Scores, s2Scores);

        double similarityScore, g1, g2, h1, h2;

        for (int i = 1, k = s2Length; i < s1Length; i++, k += s2Length) {
            h = Double.NEGATIVE_INFINITY;
            vDiagonal = s2Scores[0];
            for (int j = 1, l = k + 1; j < s2Length; j++, l++) {
                similarityScore = getSimilarityScore(s1Chars[i - 1], s2Chars[j - 1]);
                // Fill the matrices
                f = vDiagonal + similarityScore;
                g1 = s1Scores[j] + gapExtensionPenalty;
                g2 = s2Scores[j] + gapOpeningPenalty;

                if (g1 > g2) {
                    s1Scores[j] = g1;
                } else {
                    s1Scores[j] = g2;
                }

                h1 = h + gapExtensionPenalty;
                h2 = s2Scores[j - 1] + gapOpeningPenalty;
                if (h1 > h2) {
                    h = h1;
                } else {
                    h = h2;
                }

                vDiagonal = s2Scores[j];
                s2Scores[j] = Math.max(Math.max(f, s1Scores[j]), Math.max(h, 0));//maximum(f, seqOneScores[j], h, 0);

                // Set the traceback start at the current cell i, j and score
                score = Math.max(score, s2Scores[j]);

            }
        }
        return score;
    }

    private void init(int seqTwoLength, double[] seqOneScores, double[] seqTwoScores) {
        for (int jjj = 1; jjj < seqTwoLength; jjj++) {
            seqOneScores[jjj] = Double.NEGATIVE_INFINITY;
            seqTwoScores[jjj] = 0;
        }
    }

    /*
     * http://en.wikipedia.org/wiki/BLOSUM
     */
    private int[][] blosum62 =

            {{4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4},
                    {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4},
                    {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4},
                    {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4},
                    {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
                    {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4},
                    {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
                    {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4},
                    {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4},
                    {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4},
                    {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4},
                    {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4},
                    {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4},
                    {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4},
                    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
                    {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4},
                    {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4},
                    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4},
                    {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4},
                    {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4},
                    {-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4},
                    {-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
                    {0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4},
                    {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1}};

    private int getSimilarityScore(char aC, char bC) {
        int locA = charToLoc(aC);
        int locB = charToLoc(bC);
        return blosum62[locA][locB];
    }

    private int charToLoc(char a) {

        switch (a) {
            case 'A':   //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 0;
            case 'R':   //R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 1;
            case 'N':   //N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 2;
            case 'D':   //D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 3;
            case 'C':   //C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 4;
            case 'Q':   //Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 5;
            case 'E':   //E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 6;
            case 'G':   //G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 7;
            case 'H':   //H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 8;
            case 'I':   //I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 9;
            case 'L':   //L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 10;
            case 'K':   //K  M  F  P  S  T  W  Y  V  B  Z  X  *
                return 11;
            case 'M':   //M  F  P  S  T  W  Y  V  B  Z  X  *
                return 12;
            case 'F':   //F  P  S  T  W  Y  V  B  Z  X  *
                return 13;
            case 'P':   //P  S  T  W  Y  V  B  Z  X  *
                return 14;
            case 'S':   //S  T  W  Y  V  B  Z  X  *
                return 15;
            case 'T':   //S  T  W  Y  V  B  Z  X  *
                return 16;
            case 'W':   //W  Y  V  B  Z  X  *
                return 17;
            case 'Y':   //Y  V  B  Z  X  *
                return 18;
            case 'V':   //V  B  Z  X  *
                return 19;
            case 'B':   //B  Z  X  *
                return 20;
            case 'Z':   //Z  X  *
                return 21;
            case 'X':   //X  *
                return 22;
            default:
                throw new IllegalArgumentException("Not a standard amino acid!" + a);
        }

    }

}

