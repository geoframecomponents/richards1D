package richards_classes;

public class SimpleBoundaryConditionFactory {
	public BoundaryCondition createBoundaryCondition (String type) {

		BoundaryCondition boundCond = null;
		if(type.equalsIgnoreCase("Top Dirichlet") || type.equalsIgnoreCase("TopDirichlet")){
			boundCond = new TopBoundaryConditionDirichlet();
		}
		else if(type.equalsIgnoreCase("Top Neumann") || type.equalsIgnoreCase("TopNeumann")){
			boundCond = new TopBoundaryConditionNeumann();
		}
		else if(type.equalsIgnoreCase("Bottom Free Drainage") || type.equalsIgnoreCase("BottomFreeDrainage")){
			boundCond = new BottomBoundaryConditionFreeDrainage();
		}
		else if(type.equalsIgnoreCase("Bottom Dirichlet") || type.equalsIgnoreCase("BottomDirichlet")){
			boundCond = new BottomBoundaryConditionDirichlet();
		}
		return boundCond;
		}
}
