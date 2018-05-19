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

// Utility to factor out similar code dealing with common JComponent properties
interface IExtractor
{
	int value(IItem item);
}

abstract class AbstractExtractor implements IExtractor
{
	protected AbstractExtractor()
	{
	}
}

final class MinWidthExtractor extends AbstractExtractor
{
	static final IExtractor INSTANCE = new MinWidthExtractor();
	
	@Override public int value(IItem item)
	{
		return item.minimumWidth();
	}
}

final class PrefWidthExtractor implements IExtractor
{
	static final IExtractor INSTANCE = new PrefWidthExtractor();
	
	@Override public int value(IItem item)
	{
		return item.preferredWidth();
	}
}

final class PrefHeightExtractor implements IExtractor
{
	static final IExtractor INSTANCE = new PrefHeightExtractor();
	
	@Override public int value(IItem item)
	{
		return item.preferredHeight();
	}
}

final class BaselineExtractor implements IExtractor
{
	static final IExtractor INSTANCE = new BaselineExtractor();
	
	@Override public int value(IItem item)
	{
		return item.baseline();
	}
}
