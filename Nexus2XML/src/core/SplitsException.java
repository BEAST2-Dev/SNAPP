/**
 * @version $Id: SplitsException.java,v 1.4 2004-08-24 17:37:07 huson Exp $
 *
 * @author Daniel Huson
 */

package core;


/**
 * General splitstree exception
 */
public class SplitsException extends Exception {
    public SplitsException() {
        super();
    }

    public SplitsException(String str) {
        super(str);
    }
}

// EOF

