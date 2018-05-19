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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JLabel;

final class PlatformHelper
{
	private PlatformHelper()
	{
	}
	
	static int getDefaultAlignment()
	{
		switch (platform())
		{
			case MACINTOSH:
			return JLabel.TRAILING;
			
			case LINUX:
			case WINDOWS:
			case OTHER:
			default:
			return JLabel.LEADING;
		}
	}

	static List<BarRowItem> extractLeftItems(List<BarRowItem> items)
	{
		// Make sure tags order is correctly initialized
		initButtonsTagOrder();
		return extractItems(items, _leftTags);
	}
	
	static List<BarRowItem> extractCenterItems(List<BarRowItem> items)
	{
		// Make sure tags order is correctly initialized
		initButtonsTagOrder();
		return extractItems(items, _centerTags);
	}
	
	static List<BarRowItem> extractRightItems(List<BarRowItem> items)
	{
		// Make sure tags order is correctly initialized
		initButtonsTagOrder();
		return extractItems(items, _rightTags);
	}
	
	static private List<BarRowItem> extractItems(List<BarRowItem> items, String tags)
	{
		// Get the left part of the tags order for the current platform
		List<BarRowItem> extractedItems = new ArrayList<BarRowItem>(items.size());
		for (int i = 0; i < tags.length(); i++)
		{
			char code = tags.charAt(i);
			if (code == '_')
			{
				extractedItems.add(null);
			}
			else
			{
				addItems(extractedItems, items, Tag.fromCode(code));
			}
		}
		// Remove initial, last and duplicate consecutive gaps
		removeExtraGaps(extractedItems);
		return extractedItems;
	}

	static private void removeExtraGaps(List<BarRowItem> items)
	{
		// First remove duplicates
		Iterator<BarRowItem> iterator = items.iterator();
		BarRowItem previous = null;
		boolean first = true;
		while (iterator.hasNext())
		{
			BarRowItem current = iterator.next();
			if (current == null && previous == null && !first)
			{
				// Remove duplicate
				iterator.remove();
			}
			previous = current;
			first = false;
		}
		// Then remove heading and trailing gaps if any
		if (!items.isEmpty() && items.get(0) == null)
		{
			items.remove(0);
		}
		if (!items.isEmpty() && items.get(items.size() - 1) == null)
		{
			items.remove(items.size() - 1);
		}
	}
	
	static private void addItems(List<BarRowItem> target, List<BarRowItem> source, Tag tag)
	{
		boolean previousItemWasAdded = false;
		for (BarRowItem item: source)
		{
			if (item == null)
			{
				if (previousItemWasAdded)
				{
					target.add(null);
				}
				previousItemWasAdded = false;
			}
			else if (item.tag() == tag)
			{
				target.add(item);
				previousItemWasAdded = true;
			}
			else
			{
				previousItemWasAdded = false;
			}
		}
	}

	static private void initButtonsTagOrder()
	{
		if (_leftTags == null)
		{
			String tagsOrder;
			switch (platform())
			{
				case MACINTOSH:
				tagsOrder = "L_H/X/NY<>F_COA_R";
				break;
				
				case LINUX:
				tagsOrder = "L_H//XNYAC<>FO_R";
				break;

				case WINDOWS:
				case OTHER:
				default:
				tagsOrder = "L_/X/YN<>F_OCAH_R";
				break;
			}
			// Now parse into 3 strings of tags: left, center and right
			initButtonsTagOrder(tagsOrder);
		}
	}

	//CSOFF: MagicNumber
	// Used as SPI for adding new platforms, not part of "official" API
	// NB: order must be well-formed, otherwise results are unpredictable:
	// - only authorized characters
	// - no consecutive '-' duplicates
	// - exactly 2 '/' in the list
	static void initButtonsTagOrder(String order)
	{
		Pattern parser = Pattern.compile(PARSER);
		Matcher matcher = parser.matcher(order);
		if (!matcher.matches())
		{
			// This should normally never happen
			throw new IllegalArgumentException(
				"Order `" + order + "` doesn't follow required pattern");
		}
		_leftTags = group(matcher, 1);
		_centerTags = group(matcher, 2);
		_rightTags = group(matcher, 3);
	}
	//CSON: MagicNumber
	
	static private String group(Matcher matcher, int group)
	{
		String result = matcher.group(group);
		return (result == null ? "" : result);
	}
	
	static private Platform platform()
	{
		String os = System.getProperty("os.name");
		if (os.startsWith("Mac OS"))
		{
			return Platform.MACINTOSH;
		}
		else if (os.startsWith("Linux"))
		{
			return Platform.LINUX;
		}
		else if (os.startsWith("Windows"))
		{
			return Platform.WINDOWS;
		}
		else
		{
			return Platform.OTHER;
		}
	}
	
	static private enum Platform
	{
		WINDOWS,
		MACINTOSH,
		LINUX,
		OTHER
	}

	static final private String BUTTON_TAGS = "LRHYN<>FACOX_";
	static final private String TAGS_GROUP = "([" + BUTTON_TAGS + "]*)";
	static final private String PARSER = 
		"^" + TAGS_GROUP + "/" + TAGS_GROUP + "/" + TAGS_GROUP + "$";

	static private String _leftTags = null;
	static private String _centerTags = null;
	static private String _rightTags = null;
}
