/**
 * This package contains those classes necessary to integrate the energy equation. There is an abstract class and a factory
 * class that will allows to the user to choose:
 * - do not solve energy equation 
 * - solve energy equation taking into account of just pure diffusion
 * - solve energy equation taking into account of both diffusion and convection.
 * 
 * Pure diffusion is solved with an semi-implicit scheme (Picard iteration)
 * Convection is discretized  in explicit way using up wind fluxes.
 */
/**
 * @author Niccolo` Tubini
 *
 */
package Energy1DSolver;