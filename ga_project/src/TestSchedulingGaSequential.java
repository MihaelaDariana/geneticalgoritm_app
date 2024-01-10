import org.jgap.*;
import org.jgap.impl.DefaultConfiguration;
import org.jgap.impl.IntegerGene;

class TestSchedulingGASequential {

    private static final int[] TEST_DURATIONS = {15, 20, 30, 26, 22, 18, 27, 31, 24, 10};
    private static final int[] RESOURCE_REQUIREMENTS = {1, 2, 1, 3, 4, 3, 3, 2, 1, 4};

    public static void main(String[] args) throws Exception {
        Configuration conf = new DefaultConfiguration();
        conf.setPreservFittestIndividual(true);

        Gene[] genes = new Gene[TEST_DURATIONS.length];
        for (int i = 0; i < TEST_DURATIONS.length; i++) {
            genes[i] = new IntegerGene(conf, 0, TEST_DURATIONS.length - 1);
        }

        Chromosome chromosome = new Chromosome(conf, genes);
        conf.setSampleChromosome(chromosome);

        FitnessFunction fitnessFunction = new TestSchedulingFitnessFunctionSequential();
        conf.setFitnessFunction(fitnessFunction);

        conf.setPopulationSize(50); // Set the population size here

        Genotype population = Genotype.randomInitialGenotype(conf);

        for (int i = 0; i < 100; i++) {
            population.evolve();
        }

        IChromosome bestSolution = population.getFittestChromosome();

        int[] testSequence = new int[TEST_DURATIONS.length];
        int[] resourceCount = new int[RESOURCE_REQUIREMENTS.length];

        for (int i = 0; i < bestSolution.size(); i++) {
            int testIndex = (Integer) bestSolution.getGene(i).getAllele();

            int minIndex = -1;
            int minRequirement = Integer.MAX_VALUE;

            for (int j = 0; j < RESOURCE_REQUIREMENTS.length; j++) {
                if (RESOURCE_REQUIREMENTS[j] >= testIndex + 1 && RESOURCE_REQUIREMENTS[j] < minRequirement
                        && resourceCount[j] < (j + 1)) {
                    minIndex = j;
                    minRequirement = RESOURCE_REQUIREMENTS[j];
                }
            }

            if (minIndex != -1) {
                testSequence[i] = minIndex;
                resourceCount[minIndex]++;
            }
        }

        // Printing the sequence of resources assigned to tests
        for (int i = 0; i < TEST_DURATIONS.length; i++) {
            System.out.println("Test " + i + ": Duration = " + TEST_DURATIONS[i] + ", Resource Requirement = "
                    + RESOURCE_REQUIREMENTS[testSequence[i]]);
        }
    }

    static class TestSchedulingFitnessFunctionSequential extends FitnessFunction {
        @Override
        protected double evaluate(IChromosome chromosome) {
            double fitness = 0.0;
            int totalTime = 0;
            int[] resourceUsage = new int[RESOURCE_REQUIREMENTS.length];

            for (int i = 0; i < chromosome.size(); i++) {
                int testIndex = (Integer) chromosome.getGene(i).getAllele();
                totalTime += TEST_DURATIONS[testIndex];

                resourceUsage[testIndex]++;

                if (resourceUsage[testIndex] > RESOURCE_REQUIREMENTS[testIndex]) {
                    return 1.0 / (totalTime + 1); // Penalize by increasing the waiting time
                }
            }
            fitness = 1.0 / (totalTime + 1); // Return the inverse of total time

            return fitness;
        }
    }
}
