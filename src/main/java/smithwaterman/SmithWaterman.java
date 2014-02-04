

/**
 * Smith-Waterman implementation used for sequence alingments. 
 * Calculate the score of the local alignment but no alignment. This implementation uses the improvements
 * suggested by Gotoh (1989). 
 */
public class SmithWaterman {

    private double gapOpeningPenalty;
    private double gapExtensionPenalty;
    private Matrix similarityMatrx = new Blosum62()

    public SmithWaterman(double gapOpeningPenalty, double gapExtensionPenalty) {
        this.gapOpeningPenalty = gapOpeningPenalty;
        this.gapExtensionPenalty = gapExtensionPenalty;
        this.similarityMatrx = new Blosum62();

    public double run(String s1, String s2) {

        double score = 0;
        char[] seqOneChars = sequenceOne.getValue().toUpperCase().toCharArray();
        char[] seqTwoChars = sequenceTwo.getValue().toUpperCase().toCharArray();

        int seqOneLength = sequenceOne.getLength() + 1;
        int seqTwoLength = sequenceTwo.getLength() + 1;

        double f; // score of alignment x1...xi to y1...yi if xi aligns to yi
        double[] seqOneScores = new double[seqTwoLength]; // score if xi aligns to a gap after yi
        double h; // score if yi aligns to a gap after xi
        double[] seqTwoScores = new double[seqTwoLength]; // best score of alignment x1...xi to
        // y1...yi
        double vDiagonal;

        seqOneScores[0] = Double.NEGATIVE_INFINITY;
        seqTwoScores[0] = 0;

        init(seqTwoLength, seqOneScores, seqTwoScores);

        double similarityScore, g1, g2, h1, h2;

        for (int i = 1, k = seqTwoLength; i < seqOneLength; i++, k += seqTwoLength) {
            h = Double.NEGATIVE_INFINITY;
            vDiagonal = seqTwoScores[0];
            for (int j = 1, l = k + 1; j < seqTwoLength; j++, l++) {
                similarityScore = similarityMatrx.getScore(seqOneChars[i - 1], seqTwoChars[j - 1]);
                // Fill the matrices
                f = vDiagonal + similarityScore;
                g1 = seqOneScores[j] + gapExtensionPenalty;
                g2 = seqTwoScores[j] + gapOpeningPenalty;

                if (g1 > g2) {
                    seqOneScores[j] = g1;
                } else {
                    seqOneScores[j] = g2;
                }

                h1 = h + gapExtensionPenalty;
                h2 = seqTwoScores[j - 1] + gapOpeningPenalty;
                if (h1 > h2) {
                    h = h1;
                } else {
                    h = h2;
                }

                vDiagonal = seqTwoScores[j];
                seqTwoScores[j] = Math.max(Math.max(f, seqOneScores[j]), Math.max(h, 0));//maximum(f, seqOneScores[j], h, 0);

                // Set the traceback start at the current cell i, j and score
                score = Math.max(score, seqTwoScores[j]);

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

}
