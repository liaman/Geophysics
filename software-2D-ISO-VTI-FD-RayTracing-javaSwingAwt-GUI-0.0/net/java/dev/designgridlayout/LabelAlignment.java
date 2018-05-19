//  Copyright 2005-2011 Jason Aaron Osgood, Jean-Francois Poilpret
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package net.java.dev.designgridlayout;

import javax.swing.JLabel;

/**
 * Defines how components located in a sub-grid label column can be aligned inside
 * this column. 
 *
 * @see DesignGridLayout#labelAlignment(LabelAlignment)
 * @author Jean-Francois Poilpret
 */
public enum LabelAlignment
{
	/**
	 * Components in label column should be left-aligned; note that "left" means
	 * left only in left-to-right Locale orientation, otherwise it means right.
	 */
	LEFT(JLabel.LEADING),
	
	/**
	 * Components in label column should be right-aligned; note that "right" means
	 * right only in left-to-right Locale orientation, otherwise it means left.
	 */
	RIGHT(JLabel.TRAILING),

	/**
	 * Components in label column should be aligned according to guidelines for the
	 * current platform; eg on MacOS, platform alignment is right, whereas on 
	 * Windows, it is left.
	 */
	PLATFORM(PlatformHelper.getDefaultAlignment());
	
	private LabelAlignment(int align)
	{
		_align = align;
	}
	
	int alignment()
	{
		return _align;
	}
	
	final private int _align;
}
