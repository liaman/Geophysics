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

import java.awt.Component;
import java.awt.Dimension;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JComponent;
import javax.swing.JProgressBar;
import javax.swing.JSeparator;

// This helper solves an issue with baselines in Windows XP LAF on Java 5
final class BaselineHelper
{
	private BaselineHelper()
	{
	}

	static public int getBaseline(Component comp)
	{
		Dimension size = comp.getPreferredSize();
		int baseline = getBaseline(comp, size.width, size.height);
		if (baseline < 0)
		{
			boolean isCenter = false;
			// Special fix for some components with -1 baselines
			for (Class<? extends JComponent> clazz: _centerAlignedComponents)
			{
				if (clazz.isInstance(comp))
				{
					isCenter = true;
					break;
				}
			}
			if (!isCenter)
			{
				baseline = 0;
			}
		}
		return baseline;
	}

	static private int getBaseline(Component comp, int width, int height)
	{
		return comp.getBaseline(width, height);
	}

	static private final Set<Class<? extends JComponent>> _centerAlignedComponents =
		new HashSet<Class<? extends JComponent>>();
	
	static
	{
		_centerAlignedComponents.add(JSeparator.class);
		_centerAlignedComponents.add(JProgressBar.class);
	}
}
