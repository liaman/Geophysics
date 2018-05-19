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

import java.awt.Color;

import javax.swing.JComponent;
import javax.swing.JLabel;

final class MarkerHelper
{
	private MarkerHelper()
	{
	}
	
	static JComponent createMarker(int span, String tooltip)
	{
		JLabel marker = new JLabel(MARKER_LABEL);
		marker.setName(MARKER_NAME);
		marker.setHorizontalAlignment(JLabel.CENTER);
		marker.setOpaque(true);
		marker.setBackground(Color.RED);
		marker.setToolTipText(tooltip);
		return marker;
	}
	
	static boolean isMarker(JComponent component)
	{
		return		(component instanceof JLabel)
				&&	(MARKER_NAME.equals(component.getName()));
	}

	static final private String MARKER_LABEL = "spanRow()"; 
	static final private String MARKER_NAME = "DesignGridLayout.spanRow-Marker"; 
}
