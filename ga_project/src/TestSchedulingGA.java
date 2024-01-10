import org.jgap.*;
import org.jgap.impl.DefaultConfiguration;
import org.jgap.impl.IntegerGene;

public class TestSchedulingGA {

    private static final int[] TEST_DURATIONS = {15, 20, 30, 26, 22, 18, 27, 31, 24, 10};
    private static final int[] RESOURCE_REQUIREMENTS = {1, 2, 1, 3, 4, 3, 3, 2, 1, 4};
//        private static final int[] TEST_DURATIONS = {5, 7, 3, 6};
//        private static final int[] RESOURCE_REQUIREMENTS = {1, 2, 1, 3};

    public static void main(String[] args) throws Exception {
        Configuration conf = new DefaultConfiguration();
        conf.setPreservFittestIndividual(true);

        Gene[] genes = new Gene[TEST_DURATIONS.length];
        for (int i = 0; i < TEST_DURATIONS.length; i++) {
            genes[i] = new IntegerGene(conf, 0, TEST_DURATIONS.length - 1);
        }

        Chromosome chromosome = new Chromosome(conf, genes);
        conf.setSampleChromosome(chromosome);

        FitnessFunction fitnessFunction = new TestSchedulingFitnessFunction();
        conf.setFitnessFunction(fitnessFunction);

        conf.setPopulationSize(50); // Set the population size here

        Genotype population = Genotype.randomInitialGenotype(conf);

        for (int i = 0; i < 100; i++) {
            population.evolve();
        }

        IChromosome bestSolution = population.getFittestChromosome();

//        // Decode the chromosome to get the sequence of tests
//        int[] testSequence = new int[TEST_DURATIONS.length];
//        for (int i = 0; i < TEST_DURATIONS.length; i++) {
//            testSequence[i] = (Integer) bestSolution.getGene(i).getAllele();
//        }
//
//        // Print the test sequence
//        for (int i = 0; i < TEST_DURATIONS.length; i++) {
//            int testIndex = testSequence[i];
//            System.out.println("Test " + testIndex + ": Duration = " + TEST_DURATIONS[testIndex]
//                    + ", Resource Requirement = " + RESOURCE_REQUIREMENTS[testIndex]);
//        }

        int[] testSequence = new int[TEST_DURATIONS.length];
        boolean[] testPresent = new boolean[TEST_DURATIONS.length];
        for (int i = 0; i < TEST_DURATIONS.length; i++) {
            int testIndex = (Integer) bestSolution.getGene(i).getAllele();
            if (!testPresent[testIndex]) {
                testSequence[i] = testIndex;
                testPresent[testIndex] = true;
            } else {
                // Find the next available test not in the sequence
                for (int j = 0; j < TEST_DURATIONS.length; j++) {
                    if (!testPresent[j]) {
                        testSequence[i] = j;
                        testPresent[j] = true;
                        break;
                    }
                }
            }
        }

        for (int i = 0; i < TEST_DURATIONS.length; i++) {
            int testIndex = testSequence[i];
            System.out.println("Test " + testIndex + ": Duration = " + TEST_DURATIONS[testIndex]
                    + ", Resource Requirement = " + RESOURCE_REQUIREMENTS[testIndex]);
        }
    }

    static class TestSchedulingFitnessFunction extends FitnessFunction {
        @Override
//        protected double evaluate(IChromosome chromosome) {
//            double fitness = 0.0;
//            int totalTime = 0;
//            int[] resourceUsage = new int[RESOURCE_REQUIREMENTS.length];
//
//            for (int i = 0; i < chromosome.size(); i++) {
//                int testIndex = (Integer) chromosome.getGene(i).getAllele();
//                totalTime += TEST_DURATIONS[testIndex];
//
//                resourceUsage[testIndex]++;
//
//                if (resourceUsage[testIndex] > RESOURCE_REQUIREMENTS[testIndex]) {
//                    return 1.0 / (totalTime + 1); // Penalize by increasing the waiting time
//                }
//            }
//            fitness = 1.0 / (totalTime + 1); // Return the inverse of total time
//
//            return fitness;
//        }
        protected double evaluate(IChromosome chromosome) {
            int[] resourceUsage = new int[RESOURCE_REQUIREMENTS.length];
            int totalWaitingTime = 0;

            for (int i = 0; i < chromosome.size(); i++) {
                int testIndex = (Integer) chromosome.getGene(i).getAllele();
                int requiredResource = RESOURCE_REQUIREMENTS[testIndex];

                // Check if the required resource is available
                if (resourceUsage[requiredResource - 1] == 0) {
                    resourceUsage[requiredResource - 1] = testIndex + 1;
                } else {
                    // Find the next available resource that satisfies the test requirement
                    boolean resourceFound = false;
                    for (int j = requiredResource; j < RESOURCE_REQUIREMENTS.length; j++) {
                        if (resourceUsage[j] == 0) {
                            resourceUsage[j] = testIndex + 1;
                            totalWaitingTime += j - (requiredResource - 1);
                            resourceFound = true;
                            break;
                        }
                    }

                    // If no available resource found, wrap around to the beginning
                    if (!resourceFound) {
                        for (int j = 0; j < requiredResource - 1; j++) {
                            if (resourceUsage[j] == 0) {
                                resourceUsage[j] = testIndex + 1;
                                totalWaitingTime += RESOURCE_REQUIREMENTS.length - (requiredResource - 1) + j;
                                break;
                            }
                        }
                    }
                }
            }

            double fitness = 1.0 / (totalWaitingTime + 1); // Lower waiting time, higher fitness
            return fitness;
        }
    }
}
