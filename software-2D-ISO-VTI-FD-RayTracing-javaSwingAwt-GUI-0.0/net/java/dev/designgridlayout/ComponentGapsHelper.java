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

import java.awt.Container;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.LayoutStyle;
import javax.swing.LayoutStyle.ComponentPlacement;
import javax.swing.SwingConstants;

final class ComponentGapsHelper
{
	static public ComponentGapsHelper instance()
	{
		return new ComponentGapsHelper();
	}
	
	public int getHorizontalIndent()
	{
		if (_indent == -1)
		{
			JLabel label1 = new JLabel("Top label");
			JLabel label2 = new JLabel("Bottom label");
			_indent = _style.getPreferredGap(
				label1, label2, ComponentPlacement.INDENT, SwingConstants.SOUTH, null);
		}
		return _indent;
	}
	
	public int getVerticalGap(JComponent component1, JComponent component2, 
		ComponentPlacement type, Container parent)
	{
		return _style.getPreferredGap(
			component1, component2, type, SwingConstants.SOUTH, parent);
	}
	
	public int getHorizontalGap(JComponent component1, JComponent component2, 
		ComponentPlacement type, Container parent)
	{
		return _style.getPreferredGap(
			component1, component2, type, SwingConstants.EAST, parent);
	}

	public int getNorthContainerGap(JComponent component, Container parent)
	{
		return _style.getContainerGap(component, SwingConstants.NORTH, parent);
	}

	public int getSouthContainerGap(JComponent component, Container parent)
	{
		return _style.getContainerGap(component, SwingConstants.SOUTH, parent);
	}

	public int getWestContainerGap(JComponent component, Container parent)
	{
		return _style.getContainerGap(component, SwingConstants.WEST, parent);
	}

	public int getEastContainerGap(JComponent component, Container parent)
	{
		return _style.getContainerGap(component, SwingConstants.EAST, parent);
	}

	private final LayoutStyle _style = LayoutStyle.getInstance(); 
	private int _indent = -1;
}
