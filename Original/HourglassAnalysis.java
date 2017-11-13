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
		
		private HashMap<String, Double> numOfTargetPath;
		private HashMap<String, Double> numOfSourcePath;
		private HashMap<String, Double> nodePathThrough;
		private double nTotalPath;
		
		private double pathCoverageTau = 0.9;
		private HashSet<String> coreNodes;
		private HashSet<String> skipNodes;
		
		public static boolean processingFlat = false;
		private HashMap<String, Double> edgeWeights;
		
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
			addSuperSourceTarget();
			getPathStats();
			getCore();
			if (!processingFlat) {
				getFlatNetwork();
			}
		}
		
		private void addEdge(String server, String dependent) {
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
		
		private void removeEdge(String server, String dependent) {
			serves.get(server).remove(dependent);
			depends.get(dependent).remove(server);
		}
		
		private void addSuperSourceTarget() {
			for (String s: sources) {
				addEdge("superSource", s);
				if (processingFlat) {
					edgeWeights.put("superSource#" + s, 1.0);
				}
			}
			for (String t: targets) {
				addEdge(t, "superTarget");
				if (processingFlat) {
					edgeWeights.put(t + "#superTarget", 1.0);
				}
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
				addEdge(server, dependent);
			}
			scanner.close();
		}
		
		private void sourcePathsTraverse(String node) {
			if (numOfSourcePath.containsKey(node)) { // node already traversed
				return;
			}
			double nPath = 0;
			if (!skipNodes.contains(node)) {
				if (node.equals("superSource")) {
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
				if (node.equals("superTarget")) {
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
				if (numOfSourcePath.containsKey(s) && numOfTargetPath.containsKey(s)) {
					nPath = numOfSourcePath.get(s) * numOfTargetPath.get(s);
				}
				nodePathThrough.put(s, nPath);
			}
			nTotalPath = nodePathThrough.get("superSource");
		}
		
		private boolean isSuperNode(String n) {
			if (n.equals("superSource") || n.equals("superTarget")) {
				return true;
			}
			return false;
		}
			
		public void printNetworkProperties() throws Exception {
//			PrintWriter pw = new PrintWriter(new File("hourglassAnalysis.txt"));
			for (String s: nodes) {
				if (isSuperNode(s)) continue;
				System.out.println(s + "\tComplexity: " + numOfSourcePath.get(s) + "\tGenerality: " + numOfTargetPath.get(s) + "\tPath_centrality: " + nodePathThrough.get(s));
//				System.out.println(s + "\t" + nodePathThrough.get(s));
			}
			System.out.println("Total_paths: " + nTotalPath);
			System.out.println("Core_size: " + coreNodes.size() + "\t" + "Core_set: " + coreNodes);
//			pw.close();
		}
		
		private void getCore() {
			greedyTraverse(0, nTotalPath);
			// resetting everything
			skipNodes.clear();
			getPathStats();
		}
		
		private void greedyTraverse(double cumulativePathCovered, double nPath) {			
			if (!(cumulativePathCovered < nPath * pathCoverageTau)) {
				coreNodes = new HashSet(skipNodes);
				return;
			}
			double maxPathCovered = -1;
			String maxPathCoveredNode = "";
			for (String s : nodes) {
				if (isSuperNode(s)) continue;
				double numPathCovered = nodePathThrough.get(s);				
				if (numPathCovered > maxPathCovered) {
					maxPathCovered = numPathCovered;
					maxPathCoveredNode = s;
				}
			}				
			skipNodes.add(maxPathCoveredNode);
			getPathStats();
			if (!processingFlat) {
				//System.out.println(maxPathCoveredNode + "\t" + ((cumulativePathCovered + maxPathCovered) / nPath));
			}
			greedyTraverse(cumulativePathCovered + maxPathCovered, nPath);
		}
		
		private void getFlatNetwork() throws Exception {
			PrintWriter pw = new PrintWriter(new File("flat.txt"));			
			for (String s: sources) {
				removeEdge("superSource", s);
			}
			for (String s: sources) {
				addEdge("superSource", s);
				getPathStats();
				for (String r: targets) {
					if (nodePathThrough.get(r) > 0) {
						pw.println(s + "\t" + r + "\t" + numOfSourcePath.get(r));							
					}
				}
				removeEdge("superSource", s);
			}
			pw.close();
			// resetting everything
			for (String s: sources) {
				addEdge("superSource", s);
			}
			getPathStats();
		}
	}
	
	public static void main(String[] args) throws Exception {
		RPModelDAG rpModelDAG = new RPModelDAG();
		rpModelDAG.generateSimpleModel(Double.parseDouble(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]));
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
