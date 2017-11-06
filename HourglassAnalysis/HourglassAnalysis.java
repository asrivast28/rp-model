import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class HourglassAnalysis {
	private static class DependencyDAG{
		private HashSet<String> nodes;
		private HashSet<String> targets;
		private HashSet<String> sources;
		private HashMap<String, HashSet<String>> serves; 
		private HashMap<String, HashSet<String>> depends;
		
		public HashMap<String, Double> numOfTargetPath;
		public HashMap<String, Double> numOfSourcePath;
		public HashMap<String, Double> nodePathThrough;
		public double nTotalPath;
		
		public double pathCoverageTau = 0.9;
		public HashSet<String> coreNodes;
		public HashSet<String> skipNodes;
		
		public static boolean processingFlat = false;
		public HashMap<String, Double> edgeWeights;
		
		public DependencyDAG() { 
			nodes = new HashSet();
			serves = new HashMap();
			depends = new HashMap();
			targets = new HashSet();
			sources = new HashSet();
			numOfTargetPath = new HashMap();
			numOfSourcePath = new HashMap();
			nodePathThrough = new HashMap();
			skipNodes = new HashSet();
			edgeWeights = new HashMap();
		}
		
		public DependencyDAG(String dependencyGraphFilePath, String sourceFilePath, String targetFilePath) throws Exception {
			this();
			loadNetwork(dependencyGraphFilePath);
			loadSources(sourceFilePath);
			loadTargets(targetFilePath);
			getPathStats();
			getCore();
			if (!processingFlat) {
				getFlatNetwork();
			}
		}
		
		private void loadTargets(String fileName) throws Exception {
			Scanner scanner = new Scanner(new File(fileName));
			while (scanner.hasNext()) {
				targets.add(scanner.next());
			}
			scanner.close();
		}
		
		private void loadSources(String fileName) throws Exception {
			Scanner scanner = new Scanner(new File(fileName));
			while (scanner.hasNext()) {
				sources.add(scanner.next());
			}
			scanner.close();
		}
		
		private void loadNetwork(String fileName) throws Exception {
			Scanner scanner = new Scanner(new File(fileName));
			while (scanner.hasNext()) {
				String line = scanner.nextLine();
				String tokens[] = line.split("\\s+");
				String server = tokens[0];
				String dependent = tokens[1];
				if (processingFlat) {
					double weight = Double.parseDouble(tokens[2]);
					edgeWeights.put(server + "#" + dependent, weight);
				}
				nodes.add(dependent);
				nodes.add(server);
				if (serves.containsKey(server)) {
					serves.get(server).add(dependent);
				} else {
					HashSet<String> hs = new HashSet();
					hs.add(dependent);
					serves.put(server, hs);
				}
				if (depends.containsKey(dependent)) {
					depends.get(dependent).add(server);
				} else {
					HashSet<String> hs = new HashSet();
					hs.add(server);
					depends.put(dependent, hs);
				}
			}
			scanner.close();
		}
		
		private boolean isSource(String node) {		
			return sources.contains(node);
		}
		
		private boolean isTarget(String node) {
			return targets.contains(node);
		}
		
		private void sourcePathsTraverse(String node) {
			if (numOfSourcePath.containsKey(node)) { // node already traversed
				return;
			}
			double nPath = 0;
			if (!skipNodes.contains(node)) {
				if (isSource(node)) {
					++nPath;
				}
				if (depends.containsKey(node)) {
					for (String s : depends.get(node)) {
						sourcePathsTraverse(s);
						double weight = 1;
						if (processingFlat) {
							weight = edgeWeights.get(s + "#" + node);
						}
						nPath += numOfSourcePath.get(s) * weight;
					}
				}
			}
			numOfSourcePath.put(node, nPath);
		}
		
		private void targetPathsTraverse(String node) {
			if (numOfTargetPath.containsKey(node)) { // node already traversed
				return;
			}
			double nPath = 0;
			if (!skipNodes.contains(node)) {
				if (isTarget(node)) {
					++nPath;
				}
				if (serves.containsKey(node)) {
					for (String s : serves.get(node)) {
						targetPathsTraverse(s);
						double weight = 1;
						if (processingFlat) {
							weight = edgeWeights.get(node + "#" + s);
						}
						nPath += numOfTargetPath.get(s) * weight;
					}
				}
			}
			numOfTargetPath.put(node, nPath);
		}
		
		private void pathStatisticsHelper() {
			numOfSourcePath.clear();
			for (String s: nodes) {
				sourcePathsTraverse(s);
			}
			numOfTargetPath.clear();
			for (String s: nodes) {
				targetPathsTraverse(s);
			}		
		}

		private void getPathStats() {
			pathStatisticsHelper();
			nTotalPath = 0;
			for (String s : nodes) {
				double nPath = 0;
				if (numOfTargetPath.containsKey(s) && numOfSourcePath.containsKey(s)) {
					nPath = numOfTargetPath.get(s) * numOfSourcePath.get(s);
				}
				nodePathThrough.put(s, nPath);
				if (isSource(s)) {
					nTotalPath += nPath;
				}
			}
		}
			
		public void printNetworkProperties() throws Exception {
//			PrintWriter pw = new PrintWriter(new File("path_centrality.txt"));
			for (String s: nodes) {
				System.out.println(s + "\tComplexity: " + numOfSourcePath.get(s) + "\tGenerality: " + numOfTargetPath.get(s) + "\tPath_centrality: " + nodePathThrough.get(s));
			}
			System.out.println("Total_paths: " + nTotalPath);
			System.out.println("Core_size: " + coreNodes.size() + "\t" + "Core_set: " + coreNodes);
//			pw.close();
		}
		
		private void getCore() {
			greedyTraverse(0, nTotalPath);
			skipNodes.clear();
			getPathStats(); // resetting everything
		}
		
		private void greedyTraverse(double cumulativePathCovered, double nPath) {			
			if (!(cumulativePathCovered < nPath * pathCoverageTau)) {
				coreNodes = new HashSet(skipNodes);
				return;
			}
			double maxPathCovered = -1;
			String maxPathCoveredNode = "";
			for (String s : nodes) {
//				find the node with largest through path
				double numPathCovered = nodePathThrough.get(s);				
				if (numPathCovered > maxPathCovered) {
					maxPathCovered = numPathCovered;
					maxPathCoveredNode = s;
				}
			}	
			skipNodes.add(maxPathCoveredNode);
//			recompute through paths for all nodes
			getPathStats();
			greedyTraverse(cumulativePathCovered + maxPathCovered, nPath);
		}
		
		private void getFlatNetwork() throws Exception {
			PrintWriter pw = new PrintWriter(new File("flat.txt"));			
			skipNodes.clear();
			for (String s: nodes) {
				if(isSource(s)) {
					skipNodes.add(s);
				}
			}
			for (String s: nodes) {
				if (isSource(s)) {
					skipNodes.remove(s);
					getPathStats();
					for (String r: nodes) {
						if (isTarget(r) && nodePathThrough.get(r) > 0) {
							pw.println(s + "\t" + r + "\t" + nodePathThrough.get(r));							
						}
					}
					skipNodes.add(s);
				}
			}
			pw.close();
			skipNodes.clear();
			getPathStats(); // resetting everything
		}
	}
	
	public static void main(String[] args) throws Exception {
		String dependencyDAGFile = "paper_links.txt";
		String sourceFile = "paper_sources.txt";
		String targetFile = "paper_targets.txt";
		DependencyDAG dependencyDAG = new DependencyDAG(dependencyDAGFile, sourceFile, targetFile);
		System.out.println("Original_network");
		dependencyDAG.printNetworkProperties();
		String flatDAGFile = "flat.txt";
		DependencyDAG.processingFlat = true;
		DependencyDAG flatDAG = new DependencyDAG(flatDAGFile, sourceFile, targetFile);
		System.out.println("Flat_network");
		flatDAG.printNetworkProperties();
		double hScore = (1.0 - dependencyDAG.coreNodes.size() * 1.0 / flatDAG.coreNodes.size());
		System.out.println("H_score: " + hScore);
	}
}
